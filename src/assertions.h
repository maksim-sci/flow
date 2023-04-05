#define assert_eq(a,b) {if(!((a)==(b))) {fmt::print("{}:{} {}!={}: {}!={}\n",__FILE__,__LINE__,#a,#b,(a),(b));exit(-1);}}
#define assert_neq(a,b) {if((a)==(b)) {fmt::print("{}:{} {}=={}: {}!={}\n",__FILE__,__LINE__,#a,#b,(a),(b));exit(-1);}}
#define assert_tst(a) {if(!(a)) {fmt::print("{}:{} ({})==false",__FILE__,__LINE__,#a); exit(-1);}}
#define assert_simple(a) {if(!(a)) {fmt::print("assertion failed: {}:{}",__FILE__,__LINE__); exit(-1);}}