#pragma once

#include <regex>

/**
 *
 * Collection of regular expressions for the model and system readers.
 *
 */

namespace MRDEF 
{
  using namespace std;

  inline const string tok       = R"(([-a-zA-Z_0-9]+))";
  inline const string filename  = R"(([-a-zA-Z_0-9.]+))";
  inline const string num = R"(([+-]?(?:(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?|[+-]?\d+[eE][+-]?\d+)))";
  inline const string sp = R"(\s+)";
  inline const string tok_op  = R"((?:)" + sp + tok + R"()?)";  // optional token (with leading space)
  inline const string filename_op  = R"((?:)" + sp + filename + R"()?)";  // optional token (with leading space)
  inline const string ini  = R"(^)";
  inline const string end  = R"($)";
                                                                //
  inline const regex RE_EMPTY ( R"(^\s*$)" );

  inline const regex RE_NUM ( num );
  inline const regex RE_STR ( tok );

  inline const regex RE_STR_STR          ( ini + tok + sp + tok + end );
  inline const regex RE_STR_STR_STR      ( ini + tok + sp + tok + sp + tok + end);
  inline const regex RE_STR_STR_STR_STR  ( ini + tok + sp + tok + sp + tok + sp + tok + end);
  inline const regex RE_STR_STR_STROPT   ( ini + tok + sp + tok + tok_op + end);

  inline const regex RE_STR_NUM_NUM      ( ini + tok + sp + num + sp + num + end);
  inline const regex RE_STR_STR_NUM      ( ini + tok + sp + tok + sp + num + end);
  inline const regex RE_STR_NUM          ( ini + tok + sp + num );

  // Sections are prefixed with a dot
  inline const regex RE_SEC              ( ini + R"(\.)" + tok + tok_op + end ); 
  inline const regex RE_SEC_NAMEOPT      ( ini + R"(\.)" + tok + filename_op + end ); 
  inline const regex RE_SEC_TOK          ( ini + R"(\.)" + tok + sp + tok + end );    
  inline const regex RE_SEC_NAME         ( ini + R"(\.)" + tok + sp + filename + end ); 
}
