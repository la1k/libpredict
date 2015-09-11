/*
  Verify that it is possible to link with the library and that the
  header files and the library have the same version number
*/

#include <predict/predict.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
  printf("Header file version: %06i\n", PREDICT_VERSION);
  printf("  Major: %i\n", PREDICT_VERSION_MAJOR);
  printf("  Minor: %i\n", PREDICT_VERSION_MINOR);
  printf("  Patch: %i\n", PREDICT_VERSION_PATCH);
  printf("  String: %s\n", PREDICT_VERSION_STRING);
  printf("Object file version: %06i\n", predict_version());
  printf("  Major: %i\n", predict_version_major());
  printf("  Minor: %i\n", predict_version_minor());
  printf("  Patch: %i\n", predict_version_patch());
  printf("  String: %s\n", predict_version_string());

  return predict_version() != PREDICT_VERSION; // 0 == success
}
