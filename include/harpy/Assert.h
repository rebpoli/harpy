
#pragma once

#ifdef ENABLE_ASSERTS
  #define ASSERT(cond, msg_stream ) \
    do { if (!(cond)) flog << "[ASSERT] " << msg_stream; } while (0)
#else
  #define ASSERT(cond, msg)  do { } while (0)
#endif



//#include <cstdlib>
//#include <iostream>
//#include <string>
//#include <thread>
//#include <execinfo.h>   // backtrace()
//#include <unistd.h>     // STDERR_FILENO


//#include <dlfcn.h>          // dladdr
//#include <cxxabi.h>         // __cxa_demangle
//#include <cstdio>           // popen, fgets, pclose
//#include <memory>
//#include <sstream>
//#include <vector>
//#include <iomanip>

//namespace assert_detail {

//  inline std::string demangle(const char* mangled)
//  {
//    if (!mangled) return {};
//    int status = 0;
//    std::unique_ptr<char, void(*)(void*)> real(
//        abi::__cxa_demangle(mangled, nullptr, nullptr, &status), std::free);
//    return (status == 0 && real) ? std::string(real.get()) : std::string(mangled);
//  }

//  inline std::string addr2line_one(void* pc, void* so_base, const char* so_path)
//  {
//    // Convert absolute PC to offset inside the object for PIE/DSO
//    auto addr = reinterpret_cast<uintptr_t>(pc);
//    auto base = reinterpret_cast<uintptr_t>(so_base);
//    std::ostringstream os;
//    os << "addr2line -e " << (so_path ? so_path : "/proc/self/exe")
//      << " -f -C -p -i " << std::hex << "0x" << (addr - base);

//    std::string cmd = os.str();
//    FILE* pipe = ::popen(cmd.c_str(), "r");
//    if (!pipe) return {};
//    char buf[512];
//    std::string out;
//    while (fgets(buf, sizeof(buf), pipe)) out += buf;
//    ::pclose(pipe);
//    // addr2line typically returns: 'func at file:line' (already demangled due to -C)
//    return out;
//  }

//  inline void print_stacktrace(std::ostream& os)
//  {
//    void* frames[128];
//    int n = ::backtrace(frames, 128);

//    for (int i = 0; i < n; ++i) {
//      Dl_info info{};
//      ::dladdr(frames[i], &info);

//      const char* so_path = info.dli_fname ? info.dli_fname : "/proc/self/exe";
//      std::string func = demangle(info.dli_sname);

//      // Compute offset from symbol start if available
//      //        std::ptrdiff_t sym_off = 0;
//      //        if (info.dli_saddr) {
//      //            sym_off = reinterpret_cast<char*>(frames[i]) -
//      //                      reinterpret_cast<char*>(info.dli_saddr);
//      //        }

//      os << "  [" << std::setw(2) << i << "] ";
//      //           << (func.empty() ? "??" : func)
//      //           << " +0x" << std::hex << sym_off << std::dec
//      //           << "  (" << so_path << ")\n";

//      // Try to get 'func at file:line' from addr2line
//      std::string fl = addr2line_one(frames[i], info.dli_fbase, so_path);
//      if (!fl.empty()) {
//        // Trim trailing newline(s)
//        while (!fl.empty() && (fl.back() == '\n' || fl.back() == '\r')) fl.pop_back();
//        //            os << "        -> " << fl << '\n';
//        os << fl << '\n';
//      }
//    }
//  }

//  // ---------------- Failure handler ----------------
//[[noreturn]] inline void handle_assert_fail(const char* expr,
//                                            const std::string& message,
//                                            const char* file,
//                                            int line,
//                                            const char* func)
//{
//    static std::mutex m;
//    const auto tid = std::this_thread::get_id();
//    {
//        std::lock_guard<std::mutex> lock(m);
//        std::cerr << "==================== ASSERTION FAILED ====================\n"
//                  << "Expression: " << expr << '\n'
//                  << "Message   : " << message << '\n'
//                  << "File      : " << file << '\n'
//                  << "Line      : " << line << '\n'
//                  << "Function  : " << func << '\n'
//                  << "Thread    : " << tid  << '\n'
//                  << "----------------------------------------------------------\n"
//                  << "Stacktrace:\n";
//        print_stacktrace(std::cerr);
//        std::cerr << "==========================================================\n";
//        std::cerr.flush();
//    }
//    std::abort();
//}

//} // namespace




//#include <cstdarg>
//#include <cstdio>
//#include <cstdlib>
//#include <csignal>
//#include <unistd.h>
//#include <execinfo.h>

//#if defined(USE_MPI) || defined(MPI_VERSION)
//  #include <mpi.h>
//  #define ASSERT_WITH_MPI 1
//#else
//  #define ASSERT_WITH_MPI 0
//#endif

//namespace assert_min_detail {

//// Demangle C++ symbols
//inline std::string demangle(const char* name)
//{
//    if (!name) return {};
//    int status = 0;
//    std::unique_ptr<char, void(*)(void*)> real{
//        abi::__cxa_demangle(name, nullptr, nullptr, &status), std::free};
//    return (status == 0 && real) ? std::string(real.get()) : std::string(name);
//}

//// Run addr2line to get file + line number
//inline std::string addr2line_one(void* pc, void* so_base, const char* so_path)
//{
//    std::ostringstream cmd;
//    cmd << "addr2line -f -C -p -e " << (so_path ? so_path : "/proc/self/exe")
//        << " 0x" << std::hex
//        << (reinterpret_cast<uintptr_t>(pc) - reinterpret_cast<uintptr_t>(so_base));

//    FILE* pipe = ::popen(cmd.str().c_str(), "r");
//    if (!pipe) return {};

//    char buf[512];
//    std::string result;
//    while (fgets(buf, sizeof(buf), pipe)) result += buf;
//    ::pclose(pipe);

//    // Trim trailing newline(s)
//    while (!result.empty() && (result.back() == '\n' || result.back() == '\r'))
//        result.pop_back();

//    return result;
//}

//// Safer stacktrace printing using dprintf
//inline void print_stacktrace_fd(int fd = 2)
//{
//    void* frames[128];
//    int n = ::backtrace(frames, 128);

//    dprintf(fd, "Stacktrace (%d frames):\n", n);

//    for (int i = 0; i < n; ++i) {
//        Dl_info info{};
//        ::dladdr(frames[i], &info);

//        const char* so_path = info.dli_fname ? info.dli_fname : "/proc/self/exe";
//        std::string func = demangle(info.dli_sname);
//        std::string file_line = addr2line_one(frames[i], info.dli_fbase, so_path);

//        // If addr2line succeeds → print it
//        if (!file_line.empty()) {
//            dprintf(fd, "  [%02d] %s\n", i, file_line.c_str());
//        }
//        // Otherwise, fallback to function + shared object path
//        else {
//            dprintf(fd, "  [%02d] %s (%s)\n",
//                    i,
//                    func.empty() ? "??" : func.c_str(),
//                    so_path);
//        }
//    }
//}


//inline void dump_backtrace_fd(int fd)
//{
//    void* frames[64];
//    int n = ::backtrace(frames, 64);
//    dprintf(fd, "Backtrace (%d frames):\n", n);
//    ::backtrace_symbols_fd(frames, n, fd); // writes directly to fd; no iostreams
//}

//inline void fail(const char* expr,
//                 const char* file,
//                 int line,
//                 const char* func,
//                 const char* fmt, ...) noexcept
//{
//    // Header
//    dprintf(2,
//        "==================== ASSERTION FAILED ====================\n"
//        "Expression: %s\n"
//        "File      : %s\n"
//        "Line      : %d\n"
//        "Function  : %s\n",
//        expr, file, line, func);
//    print_stacktrace_fd();

//#if ASSERT_WITH_MPI
//    int inited = 0, finalized = 0;
//    MPI_Initialized(&inited);
//    MPI_Finalized(&finalized);
//    if (inited && !finalized) {
//        int rank = -1, size = -1;
//        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//        MPI_Comm_size(MPI_COMM_WORLD, &size);
//        dprintf(2, "MPI      : rank %d of %d (MPI_Abort incoming)\n", rank, size);
//    }
//#endif

//    // Message
//    if (fmt) {
//        dprintf(2, "Message  : ");
//        va_list ap; va_start(ap, fmt);
//        ::vdprintf(2, fmt, ap);
//        va_end(ap);
//        [[maybe_unused]] ssize_t ignored = ::write(2, "\n", 1);

//    }

//#ifdef ASSERT_MIN_NO_BT
//    // skip backtrace if requested
//#else
//    dprintf(2, "----------------------------------------------------------\n");
//    dump_backtrace_fd(2);
//#endif
//    dprintf(2, "==========================================================\n");

//#if ASSERT_WITH_MPI
//    if (inited && !finalized) {
//        // Try to bring the whole job down cleanly
//        MPI_Abort(MPI_COMM_WORLD, 1);
//    }
//#endif
//    // Fallbacks (in case MPI_Abort returns or MPI not initialized)
//    ::raise(SIGABRT);
//    _exit(1);
//}

//} // namespace assert_min_detail
//#ifdef ENABLE_ASSERTS
//  #define ASSERT(cond, fmt, ...) do { if (!(cond)) ::assert_min_detail::fail(#cond, __FILE__, __LINE__, __func__, fmt, ##__VA_ARGS__); } while (0)

//#else
//  #define ASSERT(cond, msg)  do { } while (0)
//#endif

