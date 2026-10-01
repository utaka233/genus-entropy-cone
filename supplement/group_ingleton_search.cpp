// Exact search on all subgroup quadruples for one distinguished Ingleton pair.
// Two pair swaps preserve the functional, so i<=j and k<=l exhaust it.
// All subgroup orders and comparisons use integers. Compile with C++17.
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdexcept>
using Mask=std::array<std::uint64_t,3>;
int size(const Mask& a){return __builtin_popcountll(a[0])+__builtin_popcountll(a[1])+__builtin_popcountll(a[2]);}
Mask meet(const Mask&a,const Mask&b){return {a[0]&b[0],a[1]&b[1],a[2]&b[2]};}
int main(int argc,char**argv){
 if(argc!=2)throw std::runtime_error("usage: group_ingleton_search input.txt");
 std::ifstream f(argv[1]);int order,n;f>>order>>n;
 if(order>192||n<1)throw std::runtime_error("invalid group size");
 std::vector<Mask> h(n);for(auto&a:h)f>>a[0]>>a[1]>>a[2];
 if(!f)throw std::runtime_error("invalid masks");
 std::vector<int> sz(n),two(n*n),three(n*n*n);
 for(int i=0;i<n;i++)sz[i]=size(h[i]);
 for(int i=0;i<n;i++)for(int j=0;j<n;j++){
   auto ij=meet(h[i],h[j]);two[i*n+j]=size(ij);
   for(int k=0;k<n;k++)three[(i*n+j)*n+k]=size(meet(ij,h[k]));
 }
 std::uint64_t checked=0,offenders=0,minLeft=0,minRight=1;
 std::array<int,4> first={-1,-1,-1,-1};
 bool allJointTrivial=true;
 for(int i=0;i<n;i++)for(int j=i;j<n;j++)for(int k=0;k<n;k++)for(int l=k;l<n;l++){
   checked++;
   std::uint64_t left=std::uint64_t(sz[i])*sz[j]*two[k*n+l]*three[(i*n+j)*n+k]*three[(i*n+j)*n+l];
   std::uint64_t right=std::uint64_t(two[i*n+j])*two[i*n+k]*two[i*n+l]*two[j*n+k]*two[j*n+l];
   if(left<right){
     offenders++;
     if(first[0]<0)first={i,j,k,l};
     if(minLeft==0||left*minRight<minLeft*right){minLeft=left;minRight=right;}
     allJointTrivial &= size(meet(meet(h[i],h[j]),meet(h[k],h[l])))==1;
     if(order==120 && left*16!=right*15)throw std::runtime_error("S5 different violation ratio");
   }
 }
 std::cout<<"order "<<order<<" subgroups "<<n<<" checked "<<checked<<" offenders "<<offenders<<"\n";
 std::cout<<"first witness subgroup indices "<<first[0]<<" "<<first[1]<<" "<<first[2]<<" "<<first[3]<<"\n";
 std::cout<<"minimum LHS/RHS "<<minLeft<<"/"<<minRight<<" all offending H1234 trivial "<<allJointTrivial<<"\n";
 if(order==168&&offenders!=0)throw std::runtime_error("PSL2(7) unexpectedly violates Ingleton");
 if(order==120&&(offenders!=60||!allJointTrivial))throw std::runtime_error("S5 expected count or joint intersection fails");
}
