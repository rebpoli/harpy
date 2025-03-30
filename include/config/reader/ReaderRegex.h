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

  inline const string tok     = R"(([-a-zA-Z_0-9]+))";
  inline const string num = R"(([+-]?(?:(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?|[+-]?\d+[eE][+-]?\d+)))";
  inline const string sp = R"(\s+)";
  inline const string tok_op  = R"((?:)" + sp + tok + R"()?)";  // optional token (with leading space)
                                                                //
  inline const regex RE_EMPTY ( R"(^\s*$)" );

  inline const regex RE_NUM ( num );
  inline const regex RE_STR ( tok );
  inline const regex RE_STR_STR ( tok + sp + tok  );
  inline const regex RE_STR_STR_STR ( tok + sp + tok + sp + tok );
  inline const regex RE_STR_STR_STR_STR ( tok + sp + tok + sp + tok + sp + tok );
  inline const regex RE_STR_STR_STROPT ( tok + sp + tok + tok_op );

  inline const regex RE_STR_NUM_NUM ( tok + sp + num + sp + num );
  inline const regex RE_STR_STR_NUM ( tok + sp + tok + sp + num );
  inline const regex RE_STR_NUM ( tok + sp + num );

  // Sections are prefixed with a dot
  inline const regex sectionRE        ( R"(\.)" + tok + tok_op ); 
  inline const regex namedSectionRE   ( R"(\.)" + tok + sp + tok );    
}
