#include "count.h"

//////////////////////////////////////////////////////////////////////
// options
//////////////////////////////////////////////////////////////////////
const char* myopts = "f:k:m:l:q:";
// -f
char* fastqfile = "-";
// -k
int Kmer_Len = 15;
// -m
int min_count = 0;
// -l
float gb_limit = 0;


//////////////////////////////////////////////////////////////////////
// constants
//////////////////////////////////////////////////////////////////////
int COUNT = 0;
int LEN = 0;
int BAD_CHAR = 0;
int PRINT_SIMPLE = 1;
const char* bintoascii = "ACGT";
int bytes_per_kmer = 44; // limit size
Mer_t Forward_Mask = 0;

//  Return the binary equivalent of  ch .
unsigned  Char_To_Binary (char ch)
{
  switch  (tolower (ch))
  {
    case  'a' : return  0;
    case  'c' : return  1;
    case  'g' : return  2;
    case  't' : return  3;
    default: BAD_CHAR++; return 0;
  };

  return  0;
}


char RC(char ch)
{
  switch(toupper(ch))
  {
    case 'A': return 'T';
    case 'T': return 'A';
    case 'C': return 'G';
    case 'G': return 'C';

    default: return 'T';
  };

  return 0;
}

char NORM(char ch)
{
  switch(toupper(ch))
  {
    case 'A': return 'A';
    case 'C': return 'C';
    case 'G': return 'G';
    case 'T': return 'T';
    default: return 'A';
  };

  return 0;
}

void InitMer(Mer_t & mer)
{
  mer = 0;
}

void MerToAscii(Mer_t mer, string & s)
{
  s.erase();
  s.resize(Kmer_Len);

  for (int i = 0; i < Kmer_Len; i++)
  {
    s[Kmer_Len-i-1] = bintoascii[mer & 0x3];
    mer >>= 2;
  }
}

//  Add  ch  to  mer  on the right, sliding one character
//  off the left end of  mer .
void  Forward_Add_Ch(Mer_t & mer, char ch)
{
  mer &= Forward_Mask;
  mer <<= 2;
  mer |= Char_To_Binary (ch);
}

//  Add the Watson-Crick complement of  ch  to  mer  on the left,
//  sliding one character off the right end of  mer .
void  Reverse_Add_Ch(Mer_t & mer, char ch)
{
  mer >>= 2;
  mer |= ((long long unsigned) (3 ^ Char_To_Binary (ch)) << (2 * Kmer_Len - 2));
}


bool Fastq_Read(FILE * fp, string & s, string & hdr, string & q)
{
  int ch;

  s.erase();
  hdr.erase();
  q.erase();

  // skip till next '@' if necessary
  while  ((ch = fgetc (fp)) != EOF && ch != '@')
    ;

  if(ch == EOF)
    return false;

  // skip spaces if any
  while  ((ch = fgetc (fp)) != EOF && ch == ' ')
    ;
  if  (ch == EOF)
    return  false;
  ungetc (ch, fp);

  // put rest of line into  hdr
  while  ((ch = fgetc (fp)) != EOF && ch != '\n')
    hdr . push_back (char (ch));

  // put everything up till next '\n' into  s
  while  ((ch = fgetc (fp)) != EOF && ch != '\n')
  {
    if  (! isspace (ch))
      s . push_back (char (ch));
  }

  // skep next line
  while  ((ch = fgetc (fp)) != EOF && ch != '\n')
    ;

  // put everything up till next '\n' into  q
  while  ((ch = fgetc (fp)) != EOF && ch != '\n')
  {
    if  (! isspace (ch))
      q . push_back (char (ch));
  }
  
  if(ch == EOF)
    return false;
  else if(ch == '@')
    ungetc (ch, fp);
  
  return  true;
}
