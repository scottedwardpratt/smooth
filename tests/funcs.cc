#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>

static int noOfFiles (void)
{
  DIR *dp;
  int i = 0;
  struct dirent *ep;
  dp = opendir ("./");

  if (dp != NULL)
  {
    while (ep = readdir (dp))
      i++;

    (void) closedir (dp);
  }
  else
    perror ("Couldn't open the directory");

  return i;
}
