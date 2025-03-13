
#include "util/Json.h"
#include "util/File.h"

#include "rapidjson/istreamwrapper.h"
#include "rapidjson/error/en.h"

#include <fstream>

using namespace rapidjson;

/**
 *
 * Acha a linha de um determinado offset em um ifstream.
 *
 */
void find_line( ifstream &ifs, uint off, uint & line, uint & loff )
{
  uint ioff=0; 
  line=1;
  loff = off;
  ifs.clear(); 
  ifs.seekg(0);
  char c; 
  while (ifs.get(c)) {
    ioff++;
    if ( ioff == off ) break;
    if ( c == '\n' ) {
      line++;
      loff = off - ioff + 1;
    }
  }
}

/**
 *
 *
 */
void read_json( RJDoc & doc, const string & fn )
{
  ifstream ifs(fn);
  if (!ifs) flog << "Arquivo JSon " << fn << " não encontrado";

  IStreamWrapper isw(ifs);
  doc.ParseStream(isw);
  if (doc.HasParseError()) {
    uint off = (unsigned)doc.GetErrorOffset();
    uint line, loff;
    find_line( ifs, off, line, loff );
    flog << "Parse error " << GetParseError_En(doc.GetParseError()) << " at line " << line << ":"<<loff<<".";
  }
}
