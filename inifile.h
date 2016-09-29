/* 
   inifile.h - a small key-value table intended to hold parameters. 

   The table can be read in from an 'ini file' with the functions kvt_read_name() and kvt_read().

   The syntax of ini files should be as follows:

    1. Basic assignment of a value to a key:

         KEY = VALUE

       Spaces around the equal sign are ignore, as are leading and
       trailing spaces in KEY and VALUE.  Neither KEY nor VALUE are
       allowed to contain spaces or hash/pound signs. Any hash/pound
       sign ('#') in the VALUE is a comment, and only the characters
       before the hash are used.

    2. A line can contain zero, one or more KEY=VALUE assignments.  If
       the line contains a hash, every thing from the hash on is
       treated as a comment and is not read in.

    3. Thus, a line starting with a hash is a comment and will be
       skipped.

    4. The exception to these comment-rule is that a line starting
       with '#include ' is not a comment, but a request to include
       another input file. The name ofthat file should follow the
       '#include ', and may optionally be surrounded by double or
       single quotes.

   Implementation in inifile.c. Test case in inifiletest.c, to be called with testa.ini as an argument. 

   Ramses van Zon, SciNet, Sept 2016
 */

#ifndef INIFILEH
#define INIFILEH

#include <stdlib.h>
#include <stdio.h>

#define MAXKEYLEN 128
#define MAXVALUELEN 3968
#define MAXTABLELEN 1024

typedef struct key_value_t {
    char key[MAXKEYLEN];
    char value[MAXVALUELEN];
} key_value_t;

typedef struct key_value_table_t {
    int         length;
    key_value_t entry[MAXTABLELEN];
} key_value_table_t;

const key_value_t* kvt_lookup_entry(const key_value_table_t* table, const char* key);
const key_value_t* kvt_begin(const key_value_table_t* table);
const key_value_t* kvt_end(const key_value_table_t* table);
const key_value_t* kvt_insert(key_value_table_t* table, const char key[MAXKEYLEN], const char value[MAXVALUELEN]);

const char* kvt_lookup(const key_value_table_t* table, const char* key);
const char* kvt_lookup_with_default(const key_value_table_t* table, const char* key, const char* def);

int kvt_lookup_int(const key_value_table_t* table, const char* key);
int kvt_lookup_int_with_default(const key_value_table_t* table, const char* key, int def);

long kvt_lookup_long(const key_value_table_t* table, const char* key);
long kvt_lookup_long_with_default(const key_value_table_t* table, const char* key, long def);

long long kvt_lookup_long_long(const key_value_table_t* table, const char* key);
long long kvt_lookup_long_long_with_default(const key_value_table_t* table, const char* key, long long def);

double kvt_lookup_double(const key_value_table_t* table, const char* key);
double kvt_lookup_double_with_default(const key_value_table_t* table, const char* key, double def);

void kvt_read_name(key_value_table_t* table, const char* filename);
void kvt_read(key_value_table_t* table, FILE* f);

#endif
