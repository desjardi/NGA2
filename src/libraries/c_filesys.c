#include <sys/stat.h>
#include <sys/types.h>
void c_mkdir(const char *dirname) {mkdir(dirname,S_IRWXU|S_IRWXG|S_IRWXO);}
int  c_isdir(const char *dirname) {
   struct stat statbuf;
   if (stat(dirname,&statbuf)) {return 0;} else {return S_ISDIR(statbuf.st_mode);}
}