//  Copyright (c) 2013 University of Pennsylvania
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included in
//  all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.

/* rnapileup
     generates pileup for RNAseq data:
      - genomic strands are treated separately
      - does not support indels 
      - does not support non-continuous alignments
*/

//  1.0 - adds another column with position-along-read data
//        (Sanger encoded just like base quals; X-33 = position along
//          read from 5' end, starting at 0)

//  1.0 - Now supports soft-clipping and discards reads with indels
//
//  1.1 - Added support for paired-end reads via the --paried option
//  	2.2.1 - fixed bug when dealing with reads with length > 127 (can now handle up to 255)
//  1.2 - Fixed handling of non-continuous (spliced) alignments. Only un-spliced, continuos aligments are used for constructing pileup.

// #define DEBUGMODE

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include "sam.h"
#include "faidx.h"
#include "bam.h"
#include <deque>
#include <cctype>

using namespace std;

///////////////////////

struct Pileup {
  // this is a basic data structure used to build a queue and report pileups
  int pos; // genomic coordinate
  char ref; // reference (genomic) nucleotide
  int nreads; // read coverage (overlap)
  string pileup; // string indicating observed read nucleotides
  string quals; // string of quality scores for observed read nucleotides
  string readpos; // string of positions with reads for observed nucleotides
  char strand; // strand indicator

  // for debugging
#ifdef DEBUGMODE
  vector<string> read_ids;
#endif
  Pileup() : pos(0), ref('N'), nreads(0), pileup(), quals(), strand('+') {}
  Pileup(int p) : pos(p), ref('N'), nreads(0), pileup(), quals(), strand('+'){ }
};

void process_queue(deque<Pileup> &q, int upto_pos, bool process_all,
		   const string &ref_id) {
  // output one-based coords!!!
  while( (!q.empty()) &&
	 (process_all || (q.front().pos < upto_pos))) {
    // output all completed positions (to the left of upto)
    // OR all queued positions if chr changed (process_all == TRUE)
    //starting with the leftmost position (queue front)
    cout << ref_id << "\t"
	 << 1+(q.front().pos) << "\t"
	 << q.front().ref << "\t"
	 << q.front().nreads << "\t"
	 << q.front().pileup << "\t"
	 << q.front().quals << "\t" 
 	<< q.front().readpos<<"\tstr"<<q.front().strand;
#ifdef DEBUGMODE
    cout << "\t";
    for (int i=0; i < q.front().read_ids.size(); ++i) {
      cout << q.front().read_ids[i] << ",";
    }
#endif
      cout << "\n";
    q.pop_front();
  }
}

/////////////////////
class DNAComplementer {
  char c[256];
public:
  DNAComplementer() { 
    for(int i=0; i<256; ++i)
      c[i] = 'n';
    c['a'] = 't'; c['c'] = 'g';  c['g'] = 'c'; c['t'] = 'a';
    c['m'] = 'k'; c['r'] = 'y';  c['k'] = 'm'; c['y'] = 'r';
    c['s'] = 's'; c['w'] = 'w';
    c['v'] = 'b'; c['h'] = 'd'; c['b'] = 'v'; c['d'] = 'h';
    for(int i='a'; i<='z'; ++i)
      c[toupper(i)] = char(toupper(c[int(i)]));
  }
  char operator () (char x) const { return c[int(x)]; }
};

/////////////////////

int main(int argc, char **argv) {
  vector<string> opts;
  vector<string> args;

  for(int i=1; i < argc; ++i) {
    string a(argv[i]);
    if (a.size() > 1 && a[0] == '-')
      opts.push_back(a);
    else
      args.push_back(a);
  }
  
  bool opts_valid = true;
  bool no_ss = false;         // not strand specific
  bool pair_ends = false;    //paired end reads
  for(int i=0; i < opts.size(); ++i) {
    if (opts[i] == "--noss")
      no_ss = true;
   else	if (opts[i] == "--paired")
		pair_ends = true;
    else
      opts_valid = false;
  }

  if(args.size() < 2 || !opts_valid) {
    cerr << "USAGE: " << argv[0] << "[options] reads.bam ref.fasta\n";
    cerr << "       --noss        Not strand-specific (convert everything to +)\n	--paired		Paired end sequencing (filter anything that is not properly paired)\n";
    return 1;
  }

  string bam_fn( args[0] );
  string fas_fn( args[1] );

  // index the fasta file by finding out where each chr starts
  ifstream file(fas_fn.c_str());
  if (!file.is_open()) {
    cerr << "Failed to open FASTA file " << fas_fn << "\n";
    return 1;
  }

  faidx_t *fai = fai_load(fas_fn.c_str());


  bamFile bam_file;
  bam_header_t *bam_hdr;

  if ((bam_file = bam_open(bam_fn.c_str(), "r")) == 0) {
    cerr << "Failed to open BAM file " << bam_fn << "\n";
    return 1;
  }

  bam_hdr = bam_header_read(bam_file);
  bam1_t *bam = bam_init1();

  string curr_ref;
  string prev_ref;
  char *ref_seq = NULL;
  int ref_len(0);
  bool changed_ref = false;

  string bases(16, 'X');
  bases[1] = 'A';
  bases[2] = 'C';
  bases[4] = 'G';
  bases[8] = 'T';
  bases[15] = 'N';

  //ReverseComplementer rev_comp;

  deque<Pileup> q; // double-ended queue;
                   // holds pileups for positions being processed

  // read BAM file
  // bam variable points at current alignment record
  while( bam_read1(bam_file, bam) > 0) {
    // bam = read alignment bam1_t
    changed_ref = false;

    // skip non-unique reads
    //int num_hits = bam_aux2i( bam_aux_get(bam, "NH") );
    //if (num_hits > 1)
    //  continue;

    // skip unmapped reads
    if ((bam->core.flag & 0x4) > 0)
      continue;
	//skip if unpaired, secondary or duplicate if pair_ends is set
	if(pair_ends&&((bam->core.flag & 0x2) == 0 || bam->core.flag & 0x100 || bam->core.flag & 0x400)){
		continue;
	}	
    // check for soft-clipped reads (clipped seq is present in SEQ field but missing from reference)
    // on both ends
    bool has_indels(false);
    bool has_skip(false); // CIGAR will have N (reference skip) in it
    int nclipstart(0); // left soft-clip size 
    int nclipend(0); // right soft-clip size

    // this DOES NOT check for SPLICED alignments e.g. with N
    // should the spliced alignments be kept or filtered out
    // check BAM_CREF_SKIP!!!
    for(int i=0; i < bam->core.n_cigar; ++i) {
      if (bam_cigar_op(bam1_cigar(bam)[i]) == BAM_CINS || 
	  bam_cigar_op(bam1_cigar(bam)[i]) == BAM_CDEL )  {
	has_indels = true; 
	break;
      }
      
      if (bam_cigar_op(bam1_cigar(bam)[i]) == BAM_CREF_SKIP)
      {
          has_skip = true; // N in CIGAR string, spliced alignment!!!
          break; // skip spliced alignments!!!
      }

      if (bam_cigar_op(bam1_cigar(bam)[i]) == BAM_CSOFT_CLIP) {
	if (i == 0)
	  nclipstart = bam_cigar_oplen(bam1_cigar(bam)[i]);
	else
	  nclipend = bam_cigar_oplen(bam1_cigar(bam)[i]);
      }
    }

    // discard reads with non-continuous (spliced or with indels) alignments
    if (has_indels || has_skip)
      continue;


    // get chr name for this bam line
    string ref( bam_hdr->target_name[bam->core.tid] );

    // load chr sequence 
    if (ref != curr_ref) {
      cerr << "Encountered new ref: " << ref << "; loading...";
      if (ref_seq) {
	free(ref_seq);
	changed_ref = true;
	prev_ref = curr_ref;
      } else {
	// initialize prev_ref to ref
	prev_ref = ref;
      }
      curr_ref = ref;
      ref_seq = fai_fetch(fai, ref.c_str(), &ref_len);
      for(int i=0; i < ref_len; ++i)
	ref_seq[i] = toupper(ref_seq[i]);
      cerr << " loaded\n";
      cerr.flush();
    }
    
    int read_pos(bam->core.pos); // 0-based leftmost coordinate
    int read_len(bam->core.l_qseq); // length of the read
    string read_qual(read_len, '#');
    string read_seq(read_len, 'N');

#ifdef DEBUGMODE
    string read_id(bam1_qname(bam));
#endif

   
    read_seq.erase(0, nclipstart);
    read_qual.erase(0, nclipstart);
    read_len -= nclipstart;
    
    // process queue: output all positions to the left of current read
    // i.e. with pos < read_pos OR if chr changed, output all positions
    // currently in the queue
    process_queue(q, read_pos, changed_ref, 
		  changed_ref ? prev_ref : curr_ref);


    // build read seq
    // NOTE: VERY IMPORTANT
    //  this works only for CONTINUOUS (unspliced) alignments!!!!!
    //  spliced alignments cannot be processed here!!
    for(int i=0; i < (read_len-nclipend); ++i) {
      // bam1_seqi(s,i) returns 4bit integer representing the base
      // 1=A, 2=C, 4=G, 8=T, 15=N
      // bam1_seq returns pointer to sequence
      read_seq[i] = bases[bam1_seqi( bam1_seq(bam), nclipstart + i)];
      
      read_qual[i] = 33+bam1_qual(bam)[nclipstart + i];
    }
    //cout << "Read: " << read_seq << "\n";

    bool rev_strand = bam1_strand(bam);
    bool rev_read = false;
    if(pair_ends){
    	if(bam->core.flag & 0x80 ){
		rev_read = true;		
	}
    } 
    deque<Pileup>::iterator q_it( q.begin() );
    for(int i=0; i < (read_len-nclipend); ++i, ++q_it) {
      int g(read_pos + i);  // genomic position

      // reported POS in bam is actually the first MATCHING base, so adjust it by the starting soft clip
      //g -= nclipstart;

      if (q_it == q.end()) {
	q.push_back(Pileup(g)); // add position (read_pos+i) to the queue
	q_it = q.end()-1; // point q to the queue end
      }
      if (q_it->pos != g) {
	cerr << "ASSERT: read pos " << g
	     << " != queue pos " << q_it->pos << "\n";
#ifdef DEBUGMODE
	cerr << read_pos << " " << read_id << "\n";
#endif
	return 1;
      }

      // cout << read_seq[i] << " vs " << ref_seq[g] << "\n";


      // skip clipped bases
      //if ( i < nclipstart || i > ((read_len-nclipend)-1) )
      //	continue;
      
      if (i == 0)
	q_it->pileup += (rev_strand^rev_read && !no_ss) ? "$" : "^~";
      else if (i == ((read_len-nclipend) - 1))
	q_it->pileup += (rev_strand^rev_read && !no_ss) ? "^~" : "$";

      // make sure we don't go past end of ref seq
      if (g >= ref_len) {
	cerr << "ASSERT: genomic pos " << g << " >= " << ref_len << "\n";
	return 1;
      }
      if(rev_strand^rev_read)
            q_it->strand = '-';
      char read_nuc = read_seq[i];
      char ref_nuc = ref_seq[g];
      if(rev_read){
     	read_nuc = DNAComplementer()(read_nuc);
	ref_nuc = DNAComplementer()(ref_nuc);
      } 
      if (ref_seq[g] == read_seq[i])
	q_it->pileup += (rev_strand^rev_read && !no_ss) ? ',' : '.';
      else
	q_it->pileup += (rev_strand^rev_read && !no_ss) ? tolower(read_nuc) : read_nuc;

      q_it->readpos += char(33 + ((rev_strand^rev_read && !no_ss) ? (read_len-(1+i)) : i)); // position string

      q_it->ref = ref_nuc; // reference (genomic) nucleotide

      q_it->quals += read_qual[i]; // quality string

#ifdef DEBUGMODE
      q_it->read_ids.push_back(read_id);
#endif

      ++(q_it->nreads); // read coverage
      
    }
    
  }

  // process all remaining positions in the queue
  process_queue(q, 0, true, curr_ref);
  return 0;
} 
