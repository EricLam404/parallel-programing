#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#define base 256
#define prime 101

void search(char* pattern, char* text);

int main(int argc, char *argv[]) {
    struct stat statbuf;
    /* Check usage */
    if (argc != 3) {
        printf("Usage: %s <pattern> <file>\n", argv[0]);
        exit(1);
    }

    /* Call stat () to fill statbuf struct */ 
    if (lstat(argv[2], &statbuf ) == -1 ) {
        printf("Could not stat file %s \n" , argv[2]);
        exit(1);
    }

    char* pattern = argv[1];
    char* file_name = argv[2];
    int file_size = statbuf.st_size;
    char text[file_size];

    FILE *file = fopen(file_name, "r");
    if (file == NULL) {
        printf("Error opening file");
        exit(1);
    }

    for(int i = 0; i < file_size; i++){
        text[i] = fgetc(file);
    }

    // printf("%s \n%s \n", pattern, file_name);
    // /* Print size of file . intmax_t is a type that holds big numbers */
    // printf ("%8jd \n", (intmax_t)statbuf.st_size);
    search(pattern, text);


    fclose(file);
    return 0;
}

void search(char* pattern, char* text){
    int pattern_length = strlen(pattern);
    int hash_value = 1;
    int i, j;
    int p = 0; /* Pattern hash value */
    int t = 0; /* Text hash value */

    // pow(base, pattern_length-1)%p so it doesn't overflow
    for(i = 0; i < pattern_length - 1; i++){
        hash_value = (hash_value * base) % prime;
    }
    for(i = 0; i < pattern_length; i++){
        p = (p * base + pattern[i]) % prime;
        t = (t * base + text[i]) % prime;
    }

    for(i = 0; i <= strlen(text) - pattern_length; i++){
        if (p == t) {
            bool matches = true;
            for (j = 0; j < pattern_length; j++) {
                if (text[i + j] != pattern[j]){
                    matches = false;
                    break;
                }
            }
            if (matches) printf("%d \n", i);
        }
        t = (base * (t - text[i] * hash_value) + text[i + pattern_length]) % prime;
        
        /* if t can be negative*/
        if (t < 0) t += prime;
    }
}
