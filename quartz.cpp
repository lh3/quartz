/* Copyright (c) 2013-2014 Y. William Yu. Released under CC0 1.0 Universal. */
/* jumpgate_sparsify.cpp
 * Sparsifies the quality vector for a read based on k-mer distance from dictionary
*/

#include "global.h"
#include "jumpgate.h"

// Populates alter_this with the correct places based on the read specified in string
void compute_alter_this(bool *alter_this, char *str, read_entry_database &red) {
	readseq a;
	int lasty;
	std::vector<readseq> mer_list;
	mer_list = encode_read_vector(str);
	lasty=-100;
	for (unsigned int y=0; y<mer_list.size(); ++y) {
		if (((alter_this[y])||(y-lasty<16))&&(!(y==(mer_list.size()-1)))) {
			// Silently do nothing
		} else {
			a = mer_list[y];
			lasty=y;
			std::vector<int> snp;
			readseq ar = rev_compl(a);
			ar = rev_compl(a);
			if (a < ar) {
				snp = red.check_hamming_neighbors(a);
			} else {
				snp = red.check_hamming_neighbors(ar);
				for (unsigned int i = 0; i < snp.size(); ++i) {
					if (snp[i]>=0)
						snp[i] = 31-snp[i];
				}
			}
			if (snp.size()==0) {
				// Silently do nothing
			} else {
				std::vector<int> locations;
				locations.resize(32);
				for (auto it = snp.begin(); it!= snp.end(); ++it) {
					if (*it >= 0)
						locations[*it]=1;
				}
				for (int x=0; x<32; ++x) {
					alter_this[x+y] = alter_this[x+y] || !(locations[x]);
				}
			}
		}
	}
}

#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "sam.h"
#include "kstring.h"
#include "kseq.h"
extern "C" {
KSEQ_DECLARE(gzFile)
}

typedef struct {
	char *name, *seq, *qual;
	bam1_t *b;
} qz_record_t;

typedef struct {
	kseq_t *fq;
	BGZF *bam;
	bam_hdr_t *hdr;
	bam1_t *b;
} qz_infile_t;

typedef struct {
	kstring_t str;
	BGZF *fp;
} qz_outfile_t;

qz_infile_t *qz_open(const char *fn, int is_bam)
{
	qz_infile_t *f;
	f = (qz_infile_t*)calloc(1, sizeof(qz_infile_t));
	if (!is_bam) {
		gzFile fp;
		fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
		f->fq = kseq_init(fp);
	} else {
		f->bam = fn && strcmp(fn, "-")? bgzf_open(fn, "r") : bgzf_dopen(fileno(stdin), "r");
		f->hdr = bam_hdr_read(f->bam);
	}
	return f;
}

void qz_close(qz_infile_t *f)
{
	if (f->bam) {
		bam_destroy1(f->b);
		bam_hdr_destroy(f->hdr);
		bgzf_close(f->bam);
	}
	if (f->fq) gzclose(f->fq->f->f);
	free(f);
}

int qz_read(qz_infile_t *f, qz_record_t *r)
{
	int ret = -1;
	if (f->fq) {
		ret = kseq_read(f->fq);
		if (ret >= 0) {
			r->name = strdup(f->fq->name.s);
			r->seq = strdup(f->fq->seq.s);
			r->qual = strdup(f->fq->qual.s);
		}
	} else if (f->bam) {
		ret = bam_read1(f->bam, f->b);
		if (ret >= 0) {
			r->b = (bam1_t*)calloc(1, sizeof(bam1_t));
			bam_copy1(r->b, f->b);
		}
	}
	return ret;
}

void qz_write_free(qz_outfile_t *f, qz_record_t *r)
{
	if (!r->b) {
		bam_write1(f->fp, r->b);
		bam_destroy1(r->b);
		free(r->seq);
	} else {
		f->str.l = 0;
		kputc('@', &f->str);
		kputs(r->name, &f->str); kputc('\n', &f->str);
		kputs(r->seq,  &f->str); kputsn("\n+\n", 3, &f->str);
		kputs(r->qual, &f->str);
		puts(f->str.s);
		free(r->name); free(r->seq); free(r->qual);
	}
}

void qz_alter(qz_record_t *r, read_entry_database &red)
{
}

extern "C" {
void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);
}

typedef struct {
	qz_infile_t *in;
	qz_outfile_t *out;
	read_entry_database *red;
	int is_bam, chunk_size, n_threads, qual;
} pipeline_t;

typedef struct {
	long n_rec;
	qz_record_t *rec;
	pipeline_t *p;
} step_t;

static void worker_for(void *_data, long i, int tid) // kt_for() callback
{
	step_t *step = (step_t*)_data;
	qz_record_t *r = &step->rec[i];
	qz_alter(r, *step->p->red);
}

static void *worker_pipeline(void *shared, int step, void *in) // kt_pipeline() callback
{
	pipeline_t *p = (pipeline_t*)shared;
	if (step == 0) { // step 0: read lines into the buffer
		step_t *s;
		s = (step_t*)calloc(1, sizeof(step_t));
		s->rec = (qz_record_t*)calloc(s->n_rec, sizeof(qz_record_t));
		s->p = p;
		for (s->n_rec = 0; s->n_rec < p->chunk_size; ++s->n_rec)
			if (qz_read(p->in, &s->rec[s->n_rec]) < 0) break;
		if (s->n_rec) return s;
	} else if (step == 1) { // step 1: reverse lines
		kt_for(p->n_threads, worker_for, in, ((step_t*)in)->n_rec);
		return in;
	} else if (step == 2) { // step 2: write the buffer to output
		step_t *s = (step_t*)in;
		while (--s->n_rec > 0)
			qz_write_free(p->out, &s->rec[s->n_rec]);
		free(s->rec); free(s);
	}
	return 0;
}

static int usage(FILE *out)
{
	fprintf(stderr, "Usage: quartz [options] <dict.prefix> <in.bam>|<in.fq.gz>\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "  -q INT    convert quality of confident base to CHAR [50]\n");
	fprintf(stderr, "  -t INT    number of threads [1]\n");
	fprintf(stderr, "  -n INT    number of records in buffer [1000000]\n");
	fprintf(stderr, "  -b        the input is a BAM file\n");
	fprintf(stderr, "  -l        enable low-memory mode\n");
	return out == stdout? 0 : 1;
}

int main(int argc, char *argv[])
{
	int64checker();
	rev_compl(0);
	subst_find(0,0);

    bool lowmem = false;
	int c;
	pipeline_t p;
	memset(&p, 0, sizeof(pipeline_t));
	p.qual = 'S'; p.n_threads = 1; p.chunk_size = 1000000;
	while ((c = getopt(argc, argv, "q:t:bln:")) >= 0) {
		if (c == 'l') lowmem = true;
		else if (c == 'q') p.qual = 33 + atoi(optarg);
		else if (c == 't') p.n_threads = atoi(optarg);
		else if (c == 'n') p.chunk_size = atoi(optarg);
	}
	if (optind + 2 > argc) return usage(stderr);

	std::string dict_filename = argv[optind];
	read_entry_database red_whole(dict_filename, lowmem);
	p.red = &red_whole;
	p.in = qz_open(argv[optind+1], p.is_bam);
	p.out = (qz_outfile_t*)calloc(1, sizeof(qz_outfile_t));
	if (p.is_bam) {
		p.out->fp = bgzf_dopen(fileno(stdout), "w");
		bam_hdr_write(p.out->fp, p.in->hdr);
	}
	kt_pipeline(p.n_threads, worker_pipeline, &p, 3);
	qz_close(p.in);
	if (p.is_bam) bgzf_close(p.out->fp);
	else free(p.out->str.s);
	return 0;
}
