#if ( defined __NSPRE_MEMPROF__ )
#include <sys/resource.h>
#include <fstream>
using namespace std;

extern "C" void 
memprof( const char location[] );

void 
memprof( const char location[] ) {
#if ( defined __NSPREMEMPROF__ ) 
#if !( defined sun4_5 ) 
    ofstream memfile("memory.out", ios::app );
    struct rusage nspre_rusage;
    getrusage( RUSAGE_SELF, &nspre_rusage );
    memfile <<" Usage in KiloBytes @ "<< location<<" : " 
            << nspre_rusage.ru_maxrss << endl ;
    memfile.flush();
    memfile.close();
#endif
#endif
}
#else 
extern "C" void memprof( const char location[] );

void memprof( const char location[] ) { return; }
#endif
   
    


        
