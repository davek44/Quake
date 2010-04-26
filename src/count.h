#ifndef COUNT_H
#define COUNT_H

#include <string>
#include <vector>

//-- Include hash_map
#ifdef __GNUC__
#if __GNUC__ < 3
  #include <hash_map.h>
  namespace Sgi { using ::hash_map; }; // inherit globals
  #define HASHMAP std
#elif __GNUC__ == 3
  #include <ext/hash_map>
  #if __GNUC_MINOR__ == 0
    namespace Sgi = std;               // GCC 3.0
    #define HASHMAP std
  #else
    namespace Sgi = ::__gnu_cxx;       // GCC 3.1 and later
    #define HASHMAP __gnu_cxx
  #endif
#elif __GNUC__ > 3
  #include <ext/hash_map>
  namespace Sgi = ::__gnu_cxx;         // GCC 4.0 and later
  #define HASHMAP __gnu_cxx
#endif
#else      // ...  there are other compilers, right?
  namespace Sgi = std;
  #define HASHMAP std
#endif

using namespace::std;
using namespace::HASHMAP;

typedef long long unsigned Mer_t;
extern Mer_t Forward_Mask;

//////////////////////////////////////////////////////////////////////
// options
//////////////////////////////////////////////////////////////////////
extern const char* myopts;
// -f
extern char * fastqfile;
// -k
extern int   Kmer_Len;
// -m
extern int min_count;
// -l
extern float gb_limit;


//////////////////////////////////////////////////////////////////////
// constants
//////////////////////////////////////////////////////////////////////
extern int COUNT;
extern int LEN;
extern int BAD_CHAR;
extern int PRINT_SIMPLE;
extern const char * bintoascii;
extern int bytes_per_kmer; // limit size

//////////////////////////////////////////////////////////////////////
// methods
//////////////////////////////////////////////////////////////////////
void InitMer(Mer_t & mer);
unsigned Char_To_Binary (char ch);
void Forward_Add_Ch (Mer_t & mer, char ch);
void Reverse_Add_Ch (Mer_t & mer, char ch);
void MerToAscii(Mer_t mer, string & s);
bool Fastq_Read(FILE * fp, string & s, string & hdr, string & q);

#endif
