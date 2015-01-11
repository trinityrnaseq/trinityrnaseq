
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <signal.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>

int lock_reg(int fd, int cmd, int type, off_t offset, int whence, off_t len)
{
  struct flock lock;
  lock.l_type = type;
  lock.l_start = offset;
  lock.l_whence = whence;
  lock.l_len = len;
  return fcntl(fd, cmd, &lock);
}

pid_t lock_test(int fd, int type, off_t offset, int whence, off_t len)
{
  struct flock lock;
  int ret;
  lock.l_type = type;
  lock.l_start = offset;
  lock.l_whence = whence;
  lock.l_len = len;
  ret = fcntl(fd, F_GETLK, &lock);
  if (ret<0) {
    printf("lock_test: fcntl F_GETLK call returned %d\n", ret);
    return ret;
  }
  if (lock.l_type==F_UNLCK)
    return 0;
  // Make sure we don't return a confusing PID (in case PID not set
  // correctly by system).
  return (lock.l_pid==0) ? 1 : lock.l_pid;
}


int main(int argc, char **argv) {

  int ret = 0;
  int fid;
  if (argc<2) {
    printf("Usage: checkLock filename\n"
	   "Checks for ability to get advisory write lock on file filename.\n"
	   "Returns 0 on success, 1 on failure.\n"
	   "--checkSetLock option: Print warning if setting lock fails.\n");
    return -1;
  }
  if (0==strcmp("--checkSetLock",argv[1])) {
    // Check setting lock
    const char filename[] = "/proc/self/exe";
    int fid = open(filename, O_RDONLY, 0);
    if (fid>=0) {
      ret = lock_reg(fid, F_SETLK, F_RDLCK, 0, SEEK_SET, 0);
      if (ret<0) {
	printf("WARNING: fcntl locks are not working.\nRunning processes are not"
	       " protected against overwrite from other computers.\n"
	       "fcntl call returned %d\n", ret);
	return 0;
      }
      lock_reg(fid, F_SETLK, F_UNLCK, 0, SEEK_SET, 0);
      close(fid);
    }
    return 0;
  } else {
    // Check whether lock is set on specified file
    fid = open(argv[1], O_RDONLY, 0);
    if (fid>0) {
      ret = lock_test(fid, F_WRLCK, 0, SEEK_SET, 0);
      if (ret==1) // special return code for unknown PID case
	printf("checkLock: File %s is locked.\n", argv[1]);
      else if (ret>0)
	printf("checkLock: File %s is locked by PID %d\n", argv[1], ret);
      close(fid);
    } else {
      if ( errno != ENOENT ) {
	printf("File %s can't be opened: %s\n", argv[1], strerror(errno) );
	ret = 1;
      }
    }
    return (ret!=0 ? 1 : 0);
  }
}
