///
/// @file inifile.h
///
/// @brief A small key-value table intended to hold parameters. 
///
/// The table can be read in from an 'ini file' with the functions
/// kvt_read_name() and kvt_read().
///
/// The syntax of ini files should be as follows:
///
/// 1. Basic assignment of a value to a key:
///
///         KEY = VALUE
///
///    Spaces around the equal sign are ignore, as are leading and
///    trailing spaces in KEY and VALUE.  Neither KEY nor VALUE are
///    allowed to contain spaces or hash/pound signs. Any hash/pound
///    sign ('#') in the VALUE is a comment, and only the characters
///    before the hash are used.
///
/// 2. A line can contain zero, one or more KEY=VALUE assignments.  If
///    the line contains a hash, every thing from the hash on is
///    treated as a comment and is not read in.
///
/// 3. Thus, a line starting with a hash is a comment and will be
///    skipped.
///
/// 4. The exception to these comment-rule is that a line starting
///    with '#include ' is not a comment, but a request to include
///    another input file. The name ofthat file should follow the
///    '#include ', and may optionally be surrounded by double or
///    single quotes.
///
/// Implementation in inifile.cpp. Test case in inifiletest.c, to be
/// called with testa.ini as an argument.
///
/// @author Ramses van Zon
/// @date 2026
///
#ifndef INIFILEH
#define INIFILEH

#include <map>
#include <string>
#include <istream>

using key_value_table_t = std::map<std::string,std::string>;
using key_value_t = key_value_table_t::value_type;

key_value_t kvt_lookup_entry(const key_value_table_t& table, const std::string& key);
key_value_table_t::const_iterator kvt_begin(const key_value_table_t& table);
key_value_table_t::const_iterator kvt_end(const key_value_table_t& table);
key_value_t kvt_insert(key_value_table_t& table, const std::string& value);

std::string kvt_lookup(const key_value_table_t& table, const std::string& key);
std::string kvt_lookup_with_default(const key_value_table_t& table, const std::string& key, const std::string& def);

int kvt_lookup_bool(const key_value_table_t& table, const std::string& key);
int kvt_lookup_bool_with_default(const key_value_table_t& table, const std::string& key, bool def);

int kvt_lookup_int(const key_value_table_t& table, const std::string& key);
int kvt_lookup_int_with_default(const key_value_table_t& table, const std::string& key, int def);

long kvt_lookup_long(const key_value_table_t& table, const std::string& key);
long kvt_lookup_long_with_default(const key_value_table_t& table, const std::string& key, long def);

long long kvt_lookup_long_long(const key_value_table_t& table, const std::string& key);
long long kvt_lookup_long_long_with_default(const key_value_table_t& table, const std::string& key, long long def);

double kvt_lookup_double(const key_value_table_t& table, const std::string& key);
double kvt_lookup_double_with_default(const key_value_table_t& table, const std::string& key, double def);

void kvt_read_name(key_value_table_t& table, const std::string& filename);
void kvt_read(key_value_table_t& table, std::istream& f);

#endif
