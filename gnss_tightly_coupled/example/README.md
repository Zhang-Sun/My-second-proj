关于数据的几点说明：
1. SPP: 需要提供brd广播星历、观测文件。
2. PPP: 需要提供brd广播星历、SP3精密星历、CLK精密卫星钟、OBS观测值、ERP地球旋转参数、TID地球潮汐参数等
3. PPK：PPP基础上多了Base的OBS观测值
4. INS组合（SOL）：SOL为SPP、PPP、PPK单独解算出来的位置（x,y,z/b,l,h）因此仅需SOL文件和IMU文件
5. INS组合 ： 原始文件的组合 无论松紧，都需要上述1.2.3对应的文件
