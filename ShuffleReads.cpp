#include "htslib/kseq.h"
#include "htslib/sam.h"
#include "htslib/hts.h"

#include <pthread.h>
#include <semaphore.h>

#include <zlib.h>
KSEQ_INIT(gzFile, gzread)
#include <iostream>
#include <fstream>

#include <vector>
#include <string>
#include <algorithm>

#include <vector>
using namespace std;
int main(int argc, char* argv[]) {

	gzFile fastaFile;
	kseq_t *ks;

	fastaFile = gzopen(argv[1], "r");
	ks = kseq_init(fastaFile);
	vector<string*> reads, names;
	long space=0;
	long iter=0;
	long n=0;
	vector<int> indices;
	while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
		string *strPtr=new string;
		strPtr->assign(ks->seq.s, ks->seq.l);
		reads.push_back(strPtr);
		space+= ks->seq.l;
		iter+= ks->seq.l;
		indices.push_back(n);
		n+=1;
		if (iter > 100000000) {
			iter=0;
			cerr << space << "\t" << n << endl;
		}
		string *namePtr = new string;
		namePtr->assign(ks->name.s, ks->name.l);
		names.push_back(namePtr);
	}
	ofstream out(argv[2]);
	std::random_shuffle(indices.begin(), indices.end());
	for (int i =0; i < indices.size(); i++ ) { 
		if (reads[indices[i]]->size() > 0) {
			out << ">" << *names[indices[i]] << endl;
			out << *reads[indices[i]] << endl;
		}
	}
		
}



