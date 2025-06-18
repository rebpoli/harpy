#pragma once

#include "base/Global.h"

#include <zlib.h>
#include <string>
#include <stdexcept>

class GzStream {
    gzFile file;
    static constexpr int BUFFER_SIZE = 4096;
    char buffer[BUFFER_SIZE];

public:
    GzStream(const string& filename) 
    {
        file = gzopen(filename.c_str(), "rb");
        if (! file) flog << "Failed to open gzip file '" << filename << "'" ;
    }
    ~GzStream() { if (file) gzclose(file); }

    /**
     *
     */
    bool getline(string& line) 
    {
        if (!gzgets(file, buffer, BUFFER_SIZE)) return false;
        line = buffer;
        if (!line.empty() && line.back() == '\n') line.pop_back();
        return true;
    }

    /**
     *
     */
    explicit operator bool() const { return file != nullptr && !gzeof(file); }
};
