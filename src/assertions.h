#define assert_eq(a,b) {if(!((a)==(b))) {throw fmt::format("{}:{} {}!={}: {}!={}\n",__FILE__,__LINE__,#a,#b,(a),(b));}}
#define assert_neq(a,b) {if((a)==(b)) {throw fmt::format("{}:{} {}=={}: {}!={}\n",__FILE__,__LINE__,#a,#b,(a),(b));}}
#define assert_tst(a) {if(!(a)) {throw fmt::format("{}:{} ({})==false\n",__FILE__,__LINE__,#a);}}
#define assert_simple(a) {if(!(a)) {throw fmt::format("assertion failed: {}:{}\n",__FILE__,__LINE__);}}