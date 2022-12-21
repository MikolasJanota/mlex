//jpms:bc
/*----------------------------------------------------------------------------*\
 * File:        fmtutils.hh
 *
 * Description: Utils for parsing DIMACS-based formats.
 *
 * Author:      jpms
 * 
 * Revision:    $Id$.
 *
 *                                     Copyright (c) 2009, Joao Marques-Silva
\*----------------------------------------------------------------------------*/
//jpms:ec

#ifndef _FMTUTILS_H
#define _FMTUTILS_H 1


//jpms:bc
/*----------------------------------------------------------------------------*\
 * Utils for DIMACS format parsing
 * (This borrows **extensively** from the MiniSAT parser)
\*----------------------------------------------------------------------------*/
//jpms:ec

#define CHUNK_LIMIT 1048576

#include <zlib.h>
#include <cstdio>
#include <string>
class StreamBuffer {
    protected:
        gzFile  in;
        char    *buf;
        int     pos;
        int     size;

        void assureLookahead() {
            if (pos >= size) {
                pos  = 0;
                size = gzread(in, buf, sizeof(buf)); } }
    public:
        explicit StreamBuffer(gzFile i) : in(i), pos(0), size(0) {
            buf = new char[CHUNK_LIMIT];
            assureLookahead();
        }

        virtual ~StreamBuffer() { delete[] buf; }

        int  operator *  () { return (pos >= size) ? EOF : buf[pos]; }
        void operator ++ () { pos++; assureLookahead(); }
};

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<class B>
static void skipWhitespace(B& in) {
    while ((*in >= 9 && *in <= 13) || *in == 32)
        ++in; }

template<class B>
static void skipTrueWhitespace(B& in) {     // not including newline
    while (*in == ' ' || *in == '\t')
        ++in; }

template<class B>
static void skipTabSpace(B& in) { while (*in == 9 || *in == 32) ++in; }

template<class B>
static void skipLine(B& in) {
    for (;;){
        if (*in == EOF) return;
        if (*in == '\n') { ++in; return; }
        if (*in == '\r') { ++in; return; }
        ++in; } }

template<class B>
static bool skipEndOfLine(B& in) {     // skip newline AND trailing comment/empty lines
    if (*in == '\n' || *in == '\r') ++in;
    else             return false;
    skipComments(in);
    return true; }

template<class B>
static bool skipText(B& in, char* text) {
    while (*text != 0){
        if (*in != *text) return false;
        ++in, ++text; }
    return true; }

template<class B>
static std::string readString(B& in) {
  std::string rstr;
  skipWhitespace(in);
  while ((*in >= '0' && *in <= '9') ||
	 (*in >= 'a' && *in <= 'z') ||
	 (*in >= 'A' && *in <= 'Z') ||
	 (*in == '_') || (*in == '.')) {
      rstr.push_back(*in);
      ++in;
  }
  return rstr;
}

template<class B>
static int parseInt(B& in) {
    int    val = 0;
    bool   neg = false;
    skipWhitespace(in);
    if      (*in == '-') neg = true, ++in;
    else if (*in == '+') ++in;
    if (*in < '0' || *in > '9') fprintf(stderr, "PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
    while (*in >= '0' && *in <= '9')
        val = val*10 + (*in - '0'),
        ++in;
    return neg ? -val : val;
}

#endif /* _FMTUTILS_H */