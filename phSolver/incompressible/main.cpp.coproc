#ifndef BUILD_SHARED_LIBS
# include "pvStaticPluginsInit.h"
#endif
extern "C" int phasta ( int ,  char* argv[] );

int 
main( int argc,   
      char* argv[] ) {
#ifndef BUILD_SHARED_LIBS
  paraview_static_plugins_init();
#endif

  phasta ( argc, argv);  
  
  return 0;
}
