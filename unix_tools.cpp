#include "unix_tools.hpp"

#include <unistd.h>
#include <cstdio>
#include <cstring>


void get_RAM() {
    //output total RAM usage
    FILE *sfile = fopen("/proc/self/status", "r");
    if (!sfile)
        return;

    char line[128];
    while (fgets(line, 128, sfile) != NULL) {
        if (strncmp(line, "VmRSS:", 6) == 0) {
            printf("%s", line);
            break;
        }
    }
    fclose(sfile);
}
