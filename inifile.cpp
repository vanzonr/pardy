///
/// @file inifile.cpp
///
/// @brief Implementation of inifile.h
///
/// @author Ramses van Zon
/// @date 2026
///
#include "inifile.h"
#include <cstring>
#include <string>
#include <fstream>

key_value_t kvt_lookup_entry(const key_value_table_t& table, const std::string& key)
{
    return {key, table.at(key)};
}

std::string kvt_lookup(const key_value_table_t& table, const std::string& key)
{
    try {
        return table.at(key);
    } catch (...) {
        fprintf(stderr, "Warning: could not find a value for '%s'.\n", key.c_str());
        throw;
    }
}

std::string kvt_lookup_with_default(const key_value_table_t& table, const std::string& key, std::string def)
{
    if (table.find(key) != table.end()) 
        return table.at(key);
    else
        return def;
}

key_value_table_t::const_iterator kvt_begin(const key_value_table_t& table)
{
    return table.cbegin();
}

key_value_table_t::const_iterator kvt_end(const key_value_table_t& table)
{
    return table.cend();
}

key_value_t kvt_insert(key_value_table_t& table,
                       const std::string& key,
                       const std::string& value)
{
    table[key] = value;
    return {key, value};
}

void kvt_read_name(key_value_table_t& table, const std::string& filename);

static void TRIMWHITESPACE(char** str)
{
    if (*str && **str && (*str)[strlen(*str)-1] == '\n')
        (*str)[strlen(*str)-1] = '\0';
    while (*str && **str && ( **str == ' ' || **str == '\t' ) )
        (*str)++;
    while (*str && **str && ( (*str)[strlen(*str)-1]==' '||(*str)[strlen(*str)-1]=='\t'))
        (*str)[strlen(*str)-1]='\0';
}

static void TRIMQUOTES(char** str) {
    if ((*str)[0]=='"' && (*str)[strlen(*str)-1]=='"') {
        (*str)[strlen(*str)-1]='\0';
        (*str)++;
    }
    if ((*str)[0]=='\'' && (*str)[strlen(*str)-1]=='\'') {
        (*str)[strlen(*str)-1]='\0';
        (*str)++;
    }
}
    
static char* process_keyval_inputline(char** key, char** value)
{
    static char empty[1] = "";
    char* equalsign = strchr(*key, '=');
    if (equalsign) {
        *equalsign = '\0';
        *value = equalsign+1;
    } else {
        *value = empty;
    }
    TRIMWHITESPACE(key);
    // treat trailing # as comments
    char* hash = strchr(*value, '#');
    if (hash) *hash='\0';
    TRIMWHITESPACE(value);
    char* space = strchr(*value, ' ');
    if (!space) space = strchr(*value, '\t');
    if (space) {
        *space = '\0';
        return space+1;
    } else {
        return NULL;
    }
}

static void kvt_read_knowingly(key_value_table_t& table, std::istream& f, const std::string& openfiles)
{    
    char empty[1] = "";
    std::string cppline;
    while (! std::getline(f, cppline).eof()) {
        char* line = const_cast<char*>(cppline.c_str());
        char* key = line;
        char* value = NULL;
        if (key) {
            TRIMWHITESPACE(&key);
            if (key[0] == '#') {
                // interpret the next word as a command followed by an argument string
                // skip pound sign and immediately following white space
                char* command = key+1;
                char* args = NULL;
                while (command && ( command[0] == ' ' || command[0] == '\t' ) )
                    command++;
                // next white space or tab is delimiter
                char* space = strchr(command, ' ');
                if (!space) space = strchr(key, '\t');
                if (space) {
                    *space = '\0';
                    args = space+1;
                } else {
                    args = empty;
                }
                TRIMWHITESPACE(&args);
                // trim surrounding quotes 
                TRIMQUOTES(&args);
                if (strcmp(command,"include")==0) {
                    if (openfiles != "" && strcmp(args,openfiles.c_str())==0)
                        fprintf(stderr, "Warning: file '%s' included recursively; skipping.\n", openfiles.c_str());
                    else
                        kvt_read_name(table, args);
                }
            } else {
                while (key && key[0]) {
                    char* nextkey = process_keyval_inputline(&key, &value);
                    kvt_insert(table, key, value);
                    key = nextkey;
                }
            }
        }
    }
}

void kvt_read_name(key_value_table_t& table, const std::string& filename)
{
    std::ifstream f;
    try {
        f.open(filename);
    } catch (...) {
        fprintf(stderr, "Warning: file '%s' could not be read; skipping.\n", filename.c_str());
        throw;
    }
    kvt_read_knowingly(table, f, filename);
    f.close();
}

void kvt_read(key_value_table_t& table, std::istream& f) {
    kvt_read_knowingly(table, f, "");
}


int kvt_lookup_bool(const key_value_table_t& table, const std::string& key)
{
    std::string value = kvt_lookup(table, key);
    for (auto& ch: value) ch = std::toupper(ch);
    if (value == "FALSE" || value == "0" || value == "NO" || value == "NONE" || value =="N")
        return false;
    else
        return true;
}

int kvt_lookup_bool_with_default(const key_value_table_t& table, const std::string& key, bool def)
{
    if (table.find(key) != table.end())
        return kvt_lookup_bool(table, key);
    else
        return def;
}

int kvt_lookup_int(const key_value_table_t& table, const std::string& key)
{
    return std::stoi(kvt_lookup(table, key));
}

int kvt_lookup_int_with_default(const key_value_table_t& table, const std::string& key, int def)
{
    if (table.find(key) != table.end())
        return kvt_lookup_int(table, key);
    else
        return def;
}

long kvt_lookup_long(const key_value_table_t& table, const std::string& key)
{
    return std::stol(kvt_lookup(table, key));
}

long kvt_lookup_long_with_default(const key_value_table_t& table, const std::string& key, long def)
{
    if (table.find(key) != table.end())
        return kvt_lookup_long(table, key);
    else
        return def;
}

long long kvt_lookup_long_long(const key_value_table_t& table, const std::string& key)
{
    return std::stoll(kvt_lookup(table, key));
}

long long kvt_lookup_long_long_with_default(const key_value_table_t& table, const std::string& key, long long def)
{
    if (table.find(key) != table.end())
        return kvt_lookup_long_long(table, key);
    else
        return def;
}

double kvt_lookup_double(const key_value_table_t& table, const std::string& key)
{
    return std::stod(kvt_lookup(table, key));
}

double kvt_lookup_double_with_default(const key_value_table_t& table, const std::string& key, double def)
{
    if (table.find(key) != table.end())
        return kvt_lookup_double(table, key);
    else
        return def;
}
