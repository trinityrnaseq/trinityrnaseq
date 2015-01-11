/*
Copyright (C) 2012 Francesco Strozzi <francesco.strozzi@gmail.com>

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

#define _POSIX_C_SOURCE 1
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "kseq.h"

KSEQ_INIT(int, read)

int print_seq(int ilmn_trinity, char *append, int to_fasta, char *s[]) {
	if (ilmn_trinity && (s[0][strlen(s[0])-2] != '/') && s[3] != NULL) {
		printf(">%s/%c\n", s[0], s[3][0]);
	}
	else if (append == NULL) {
		printf(">%s\n", s[0]);
	}
	else {
		printf(">%s%s\n", s[0], append);		
	}
	printf("%s\n", s[1]);
	if (!to_fasta && s[2] != NULL && !ilmn_trinity) {
		printf("+\n%s\n",s[2]);
	}
	return 0;
}


int process_input(FILE* stream, int rev_comp, char *string, int to_fa, int ilmn_trin) {

	kseq_t *seq;
	seq = kseq_init(fileno(stream));
	int count = 0;
	int res = 0;
	if (rev_comp) {
		while ((res = kseq_read(seq)) >= 0) {
			
			char *sequence_to_print[4] = {[2] = NULL};
	
			char quality[seq->qual.l];
			if (seq->qual.s) {
				for (int i = 0; i < seq->qual.l; ++i)
				{
					quality[i] = *(seq->qual.s + seq->qual.l -1 -i);
				}
				quality[seq->qual.l] = '\0';
				sequence_to_print[2] = quality;
			}
			
			char rev_seq[seq->seq.l];
			for(int i = 0; i < seq->seq.l; ++i) {
				if (*(seq->seq.s + seq->seq.l-1 - i) == 'A' || *(seq->seq.s + seq->seq.l-1 - i) == 'a') rev_seq[i] = 'T';
				else if (*(seq->seq.s + seq->seq.l-1 - i) == 'C' || *(seq->seq.s + seq->seq.l-1 - i) == 'c') rev_seq[i] = 'G';
				else if (*(seq->seq.s + seq->seq.l-1 - i) == 'T' || *(seq->seq.s + seq->seq.l-1 - i) == 't') rev_seq[i] = 'A';
				else if (*(seq->seq.s + seq->seq.l-1 - i) == 'G' || *(seq->seq.s + seq->seq.l-1 - i) == 'g') rev_seq[i] = 'C';
				else if (*(seq->seq.s + seq->seq.l-1 - i) == 'N' || *(seq->seq.s + seq->seq.l-1 - i) == 'n') rev_seq[i] = 'N';
				else if (*(seq->seq.s + seq->seq.l-1 - i) == 'U' || *(seq->seq.s + seq->seq.l-1 - i) == 'u') rev_seq[i] = 'A';
			}
			rev_seq[seq->seq.l] = '\0';
			sequence_to_print[0] = seq->name.s;
			sequence_to_print[1] = rev_seq;
			
			if (seq->comment.s) {
				sequence_to_print[3] = seq->comment.s;
			}
			else {
				sequence_to_print[3] = NULL;
			}
			
			print_seq(ilmn_trin, string, to_fa, sequence_to_print );
			count++;
			}
		

	}
	else {
		while ((res = kseq_read(seq)) >= 0) {
			char *sequence_to_print[4];
			sequence_to_print[0] = seq->name.s;
			sequence_to_print[1] = seq->seq.s;
			for (int i=0;i < strlen(sequence_to_print[1]); i++) {
    			sequence_to_print[1][i] = toupper(sequence_to_print[1][i]);
			}
			if (seq->qual.s) sequence_to_print[2] = seq->qual.s;
			if (seq->comment.s) {
				sequence_to_print[3] = seq->comment.s;
			}
			else {
				sequence_to_print[3] = NULL;
			}
			print_seq(ilmn_trin, string, to_fa, sequence_to_print);
			count++;
		}
	}

	kseq_destroy(seq);
	fclose(stream);
	fprintf(stderr,"Sequences parsed: %d\n",count);
	if (res == -1) {
		exit(0);
	}
	else {
		fprintf (stderr,"fastool: parsing error, truncated sequence and/or quality detected!\n");
		exit(1);
	}
	
}

void print_help(char *command_line) {
		printf("Usage: %s (--rev) (--append [string_to_append_to_header]) (--to-fasta) (--illumina-trinity) sequences_1.fastq/a sequences_2.fastq/a ... \n",command_line);
}

int main(int argc, char *argv[])
{

	int reverse_complement = 0;
	char *string_to_append = NULL;
	int to_fasta = 0;
	int read_from_file = 0;
	int illumina_trinity = 0;

	if(argc == 1) {
		print_help(argv[0]);
		exit(1);
	}

	for(int i = 1; i < argc; ++i)
	{
		if(strcmp(argv[i],"--rev") == 0) reverse_complement = 1;
		else if (strcmp(argv[1],"-h") == 0) {
			print_help(argv[0]);
			exit(1);
		}
		else if(strcmp(argv[i],"--to-fasta") == 0) to_fasta = 1;
		else if((strcmp(argv[i],"--append") == 0)) {
			if (i+1 == argc) {
				printf("String to append is missing!\n");
				exit(1);
			}
			else {
				string_to_append = argv[i+1];
				i++;
			}
		}
		else if((strcmp(argv[i],"--illumina-trinity") == 0)) illumina_trinity = 1;
		else {
			if (illumina_trinity && string_to_append) {
				printf("You are using both --append and --illumina-trinity options. You can only provide one or the other.\n");
				exit(1);
			}
			read_from_file = 1;
			FILE* fp;
			if (!(fp = fopen(argv[i],"r"))) {
				printf("No %s file found!\n", argv[i]);
				exit(1);
			}
			process_input(fp,reverse_complement, string_to_append, to_fasta, illumina_trinity);
		}
	}
	if (!read_from_file)
	{
		FILE* fp_stdin;
		fp_stdin = fdopen(fileno(stdin), "rb");
		process_input(fp_stdin, reverse_complement, string_to_append, to_fasta, illumina_trinity);
	}

}



