#include "inifile.h"
#include <string.h>
#include <assert.h>

const key_value_t* kvt_lookup_entry(const key_value_table_t* table, const char* key)
{
    for (int i=0;i<table->length;i++)
        if (strcmp(key, table->entry[i].key)==0)
            return table->entry + i;
    return NULL;
}

const char* kvt_lookup(const key_value_table_t* table, const char* key)
{
    const key_value_t* entry = kvt_lookup_entry(table, key);
    if (entry != NULL)
        return entry->value;
    else {
        fprintf(stderr, "Warning: could not find a value for '%s'.\n", key);
        return NULL;
    }
}

const char* kvt_lookup_with_default(const key_value_table_t* table, const char* key, const char* def)
{
    const key_value_t* entry = kvt_lookup_entry(table, key);
    if (entry != NULL)
        return entry->value;
    else
        return def;
}

const key_value_t* kvt_begin(const key_value_table_t* table)
{
    return table->entry;
}

const key_value_t* kvt_end(const key_value_table_t* table)
{
    return table->entry + table->length;
}

const key_value_t* kvt_insert(key_value_table_t* table,
                              const char key[MAXKEYLEN],
                              const char value[MAXVALUELEN])
{
    key_value_t* existing_entry = (key_value_t*)kvt_lookup_entry(table, key);
    if (existing_entry==NULL) {
        int len = table->length;
        if (len < MAXTABLELEN) {
            strncpy(table->entry[len].key, key, MAXKEYLEN);
            strncpy(table->entry[len].value, value, MAXVALUELEN);
            table->length ++;
            return table->entry + len;
        } else {
            return NULL;
        }
    } else {        
        strncpy(existing_entry->value, value, MAXVALUELEN);
        return existing_entry;
    }
}

void kvt_read_name(key_value_table_t* table, const char* filename);

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
    /* treat trailing # as comments */
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

static void kvt_read_knowingly(key_value_table_t* table, FILE* f, const char* openfiles)
{    
    char*  line = NULL;
    size_t len = 0;
    char empty[1] = "";
    while (getline(&line, &len, f) != -1) {
        char* key = line;
        char* value = NULL;
        if (key) {
            TRIMWHITESPACE(&key);
            if (key[0] == '#') {
                /* interpret the next word as a command followed by an argument string */
                /* skip pound sign and immediately following white space */
                char* command = key+1;
                char* args = NULL;
                while (command && ( command[0] == ' ' || command[0] == '\t' ) )
                    command++;
                /* next white space or tab is delimiter */
                char* space = strchr(command, ' ');
                if (!space) space = strchr(key, '\t');
                if (space) {
                    *space = '\0';
                    args = space+1;
                } else {
                    args = empty;
                }
                TRIMWHITESPACE(&args);
                /* trim surrounding quotes */
                TRIMQUOTES(&args);
                if (strcmp(command,"include")==0) {
                    if (openfiles && strcmp(args,openfiles)==0)
                        fprintf(stderr, "Warning: file '%s' included recursively; skipping.\n", openfiles);
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
    free(line);
}

void kvt_read_name(key_value_table_t* table, const char* filename)
{
    FILE* f  = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "Warning: file '%s' could not be read; skipping.\n", filename);
    } else {
        kvt_read_knowingly(table, f, filename);
        fclose(f);
    }
}

void kvt_read(key_value_table_t* table, FILE* f) {
    kvt_read_knowingly(table, f, NULL);
}


int kvt_lookup_int(const key_value_table_t* table, const char* key)
{
    const char* value = kvt_lookup(table, key);
    if (value) {
        int result;
        sscanf(value, "%d", &result);
        return result;
    } else {
        exit(1);
    }
}

int kvt_lookup_int_with_default(const key_value_table_t* table, const char* key, int def)
{
    const char* value = kvt_lookup(table, key);
    if (value) {
        int result;
        sscanf(value, "%d", &result);
        return result;
    } else {
        return def;
    }
}

long kvt_lookup_long(const key_value_table_t* table, const char* key)
{
    const char* value = kvt_lookup(table, key);
    if (value) {
        long result;
        sscanf(value, "%ld", &result);
        return result;
    } else {
        exit(1);
    }
}

long kvt_lookup_long_with_default(const key_value_table_t* table, const char* key, long def)
{
    const char* value = kvt_lookup(table, key);
    if (value) {
        long result;
        sscanf(value, "%ld", &result);
        return result;
    } else {
        return def;
    }
}

long long kvt_lookup_long_long(const key_value_table_t* table, const char* key)
{
    const char* value = kvt_lookup(table, key);
    if (value) {
        long long result;
        sscanf(value, "%lld", &result);
        return result;
    } else {
        exit(1);
    }
}

long long kvt_lookup_long_long_with_default(const key_value_table_t* table, const char* key, long long def)
{
    const char* value = kvt_lookup(table, key);
    if (value) {
        long long result;
        sscanf(value, "%lld", &result);
        return result;
    } else {
        return def;
    }
}

double kvt_lookup_double(const key_value_table_t* table, const char* key)
{
    const char* value = kvt_lookup(table, key);
    if (value) {
        double result;
        sscanf(value, "%lf", &result);
        return result;
    } else {
        exit(1);
    }
}

double kvt_lookup_double_with_default(const key_value_table_t* table, const char* key, double def)
{
    const char* value = kvt_lookup(table, key);
    if (value) {
        double result;
        sscanf(value, "%lf", &result);
        return result;
    } else {
        return def;
    }
}
