cd .../tutorial
cmake -Dhelib_DIR=/usr/local/include/helib_pack/share/cmake/helib -DCMAKE_BUILD_TYPE=Debug .
make -j32

cd .../tutorial
cmake -Dhelib_DIR=/usr/local/include/helib_pack/share/cmake/helib -DCMAKE_BUILD_TYPE=Debug -DENABLE_GCOV=1 .
make -j32


#perf && Flame Graph
1、git将Flame Graph clone下来：git clone https://github.com/brendangregg/FlameGraph.git

2、sudo perf record -e cpu-clock --call-graph dwarf ./01_ckks_basics
sudo perf record -e cpu-clock -g ./01_ckks_basics
    sudo perf record -a -g ./01_ckks_basics
    结束执行后，在当前目录下会生成采样数据perf.data.

3、用perf script工具对perf.data进行解析
    perf script -i perf.data &> perf.unfold

4、将perf.unfold中的符号进行折叠：
    ./stackcollapse-perf.pl perf.unfold &> perf.folded

5、最后生成svg图：
    ./flamegraph.pl perf.folded > perf.svg



