#pragma once

#include "base/Global.h"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>


class CsvFile;

inline static CsvFile& endrow(CsvFile& file);
inline static CsvFile& flush(CsvFile& file);

class CsvFile
{
    std::ofstream fs_;
    bool is_first_;
    const std::string separator_;
    const std::string escape_seq_;
    const std::string special_chars_;
public:
    CsvFile(const std::string filename, const std::string separator = "\t", bool append = 1)
        : fs_()
        , is_first_(true)
        , separator_(separator)
        , escape_seq_("\"")
        , special_chars_("\"")
    {
        fs_.exceptions(std::ios::failbit | std::ios::badbit);
        if ( append )
          fs_.open(filename, std::ofstream::app);
        else
          fs_.open(filename);
    }

    ~CsvFile() { flush(); fs_.close(); }
    void flush() { fs_.flush(); }

    void endrow()
    {
        fs_ << std::endl;
        is_first_ = true;
    }

    CsvFile& operator << ( CsvFile& (* val)(CsvFile&))
    {
        return val(*this);
    }

    CsvFile& operator << (const char * val)
    {
        return write(escape(val));
    }

    CsvFile& operator << (const std::string & val)
    {
        return write(escape(val));
    }

    template<typename T>
    CsvFile& operator << (const T& val)
    {
        return write(val);
    }

protected:
    template<typename T>
    CsvFile& write (const T& val)
    {
        if (!is_first_)
        {
            fs_ << separator_;
        }
        else
        {
            is_first_ = false;
        }
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

class CsvFile1 
{
    std::ofstream fs_;
    bool is_first_;
    const std::string separator_;
    const std::string escape_seq_;
    const std::string special_chars_;
public:
    CsvFile1(const std::string filename, const std::string separator = "\t", bool append = 1)
        : fs_()
        , is_first_(true)
        , separator_(separator)
        , escape_seq_("\"")
        , special_chars_("\"")
    {
      if ( ! SINGLEPROC ) return;

      fs_.exceptions(std::ios::failbit | std::ios::badbit);
      if ( append )
        fs_.open(filename, std::ofstream::app);
      else
        fs_.open(filename);
    }

    ~CsvFile1() { if ( ! SINGLEPROC ) return; flush(); fs_.close(); }
    void flush() { fs_.flush(); }

    void endrow()
    {
        if ( ! SINGLEPROC ) return;
        fs_ << std::endl;
        is_first_ = true;
    }

    CsvFile1& operator << ( CsvFile1& (* val)(CsvFile1&))
    {
        return val(*this);
    }

    CsvFile1& operator << (const char * val)
    {
        return write(escape(val));
    }

    CsvFile1& operator << (const std::string & val)
    {
        return write(escape(val));
    }

    template<typename T>
    CsvFile1& operator << (const T& val)
    {
        return write(val);
    }

protected:
    template<typename T>
    CsvFile1& write (const T& val)
    {
        if ( ! SINGLEPROC ) return *this;
        if (!is_first_)
        {
            fs_ << separator_;
        }
        else
        {
            is_first_ = false;
        }
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


inline static CsvFile& endrow(CsvFile& file)
{
    file.endrow();
    return file;
}

inline static CsvFile& flush(CsvFile& file)
{
    file.flush();
    return file;
}

inline static CsvFile1& endrow(CsvFile1& file)
{
    file.endrow();
    return file;
}

inline static CsvFile1& flush(CsvFile1& file)
{
    file.flush();
    return file;
}
