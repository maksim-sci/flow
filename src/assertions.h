#pragma once

inline void debug_fallout()
{
    #ifdef DEBUG_MY_VERSION
    __debugbreak();
    #endif
}



#define assert_eq(a,b) {if(!((a)==(b))) {fmt::print("{}:{} {}!={}: {}!={}\n",__FILE__,__LINE__,#a,#b,(a),(b));debug_fallout();}}
#define assert_neq(a,b) {if((a)==(b)) {fmt::print("{}:{} {}=={}: {}!={}\n",__FILE__,__LINE__,#a,#b,(a),(b));debug_fallout();}}
#define assert_tst(a) {if(!(a)) {fmt::print("{}:{} ({})==false\n",__FILE__,__LINE__,#a);debug_fallout();}}
#define assert_simple(a) {if(!(a)) {fmt::print("assertion failed: {}:{}\n",__FILE__,__LINE__);debug_fallout();}}