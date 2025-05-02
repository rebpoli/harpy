#pragma once

#include "base/Global.h"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

// Formatting classes
struct CSVSci { CSVSci(double v_) : val(v_) {} ; double val; };

/**
 *
 * Only on the first processor (rank=0)
 *
 */
class CsvFile1 
{
  std::ofstream fs_;
  bool is_first_;
  const std::string separator_;
  const std::string escape_seq_;
  const std::string special_chars_;
  public:
  CsvFile1(const std::string filename, const std::string separator = "\t", bool append = 1)
    : fs_() , is_first_(true) , separator_(separator) , escape_seq_("\"") , special_chars_("\"")
  {
    if ( ! SINGLEPROC ) return;

    fs_.exceptions(std::ios::failbit | std::ios::badbit);
    if ( append )  fs_.open(filename, std::ofstream::app);
    else           fs_.open(filename);

    saveState();
  }

  ~CsvFile1() { if ( ! SINGLEPROC ) return; flush(); fs_.close(); }
  void flush() { fs_.flush(); }

  // Save/restore the current stream state
  std::ios::fmtflags savedFlags;
  std::streamsize savedPrecision;
  char savedFill;
  void saveState() { savedFlags = fs_.flags(); savedPrecision = fs_.precision(); savedFill = fs_.fill(); }
  void restoreState() { fs_.flags(savedFlags); fs_.precision(savedPrecision); fs_.fill(savedFill); }

  /**   **/
  void endrow() { if ( ! SINGLEPROC ) return; fs_ << std::endl; is_first_ = true; }
  CsvFile1& operator << ( CsvFile1& (* val)(CsvFile1&)) { return val(*this); }
  CsvFile1& operator << (const char * val) { return write(escape(val)); }
  CsvFile1& operator << (const std::string & val) { return write(escape(val)); }

  /// Formatting classes to control the stream state
  CsvFile1 & operator<<( const CSVSci & obj )
  { fs_ << scientific; write(obj.val); restoreState(); return *this; }

  /// All other types
  template<typename T> CsvFile1 & operator<< (const T& val) { return write(val); }


protected:
  CsvFile1& write (const double & val)
  {
    if ( ! SINGLEPROC ) return *this;

    if (!is_first_) { fs_ << separator_; }
    else { is_first_ = false; }

    fs_ << val;

    return *this;
  }
  template<typename T>
    CsvFile1& write (const T& val)
    {
      if ( ! SINGLEPROC ) return *this;

      if (!is_first_) fs_ << separator_;
      else is_first_ = false;

      fs_ << val;
      return *this;
    }

  std::string escape(const std::string & val)
  {
    std::ostringstream result;
    result << '"';
    std::string::size_type to, from = 0u, len = val.length();
    while (from < len &&
        std::string::npos != (to = val.find_first_of(special_chars_, from)))
    {
      result << val.substr(from, to - from) << escape_seq_ << val[to];
      from = to + 1;
    }
    result << val.substr(from) << '"';
    return result.str();
  }
};


inline static CsvFile1& endrow(CsvFile1& file) { file.endrow(); return file; }
inline static CsvFile1& flush(CsvFile1& file) { file.flush(); return file; }



/**   **/
/**
 *
 * All processors
 *
 */
class CsvFile
{
  std::ofstream fs_;
  bool is_first_;
  const std::string separator_;
  const std::string escape_seq_;
  const std::string special_chars_;
  public:
  CsvFile(const std::string filename, const std::string separator = "\t", bool append = 1)
    : fs_() , is_first_(true) , separator_(separator) , escape_seq_("\"") , special_chars_("\"")
  {
    fs_.exceptions(std::ios::failbit | std::ios::badbit);
    if ( append )  fs_.open(filename, std::ofstream::app);
    else           fs_.open(filename);

    saveState();
  }

  ~CsvFile() { flush(); fs_.close(); }
  void flush() { fs_.flush(); }

  // Save/restore the current stream state
  std::ios::fmtflags savedFlags;
  std::streamsize savedPrecision;
  char savedFill;
  void saveState() { savedFlags = fs_.flags(); savedPrecision = fs_.precision(); savedFill = fs_.fill(); }
  void restoreState() { fs_.flags(savedFlags); fs_.precision(savedPrecision); fs_.fill(savedFill); }

  /**   **/
  void endrow() { fs_ << std::endl; is_first_ = true; }
  CsvFile& operator << ( CsvFile& (* val)(CsvFile&)) { return val(*this); }
  CsvFile& operator << (const char * val) { return write(escape(val)); }
  CsvFile& operator << (const std::string & val) { return write(escape(val)); }

  /// Formatting classes to control the stream state
  CsvFile & operator<<( const CSVSci & obj )
  { fs_ << scientific; write(obj.val); restoreState(); return *this; }

  /// All other types
  template<typename T> CsvFile & operator<< (const T& val) { return write(val); }


protected:
  CsvFile& write (const double & val)
  {
    if (!is_first_) { fs_ << separator_; }
    else { is_first_ = false; }

    fs_ << val;

    return *this;
  }
  template<typename T>
    CsvFile& write (const T& val)
    {
      if (!is_first_) fs_ << separator_;
      else is_first_ = false;

      fs_ << val;
      return *this;
    }

  std::string escape(const std::string & val)
  {
    std::ostringstream result;
    result << '"';
    std::string::size_type to, from = 0u, len = val.length();
    while (from < len &&
        std::string::npos != (to = val.find_first_of(special_chars_, from)))
    {
      result << val.substr(from, to - from) << escape_seq_ << val[to];
      from = to + 1;
    }
    result << val.substr(from) << '"';
    return result.str();
  }
};

inline static CsvFile & endrow(CsvFile & file) { file.endrow(); return file; }
inline static CsvFile & flush(CsvFile & file) { file.flush(); return file; }
