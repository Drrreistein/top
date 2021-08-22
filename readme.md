## Code Usage

Requirements:`numpy` and `matplotlib` are required for executing the python script.

```shell
python3 main.py
```



## problem definition

a full connected truss structure **(3 meter * 2 meter)**: left side of the structure are fixed, its right side is loaded by **50kN**. 

| material properties  |       | unit   |
| -------------------- | ----- | ------ |
| cross section        | 200   | $mm^2$ |
| module of elasticity | 29000 | $E$    |

**object**: minimize the material usage of the structure

**constraints**:

1.  maximal nodal displacement no bigger than 0.1m
2. maximal stress of each linkage don't exceed 300 MPa

<img src="./docs/image-20210822103205089.png" alt="image-20210822103205089" style="zoom:50%;" />



## methods

<img src="./docs/topology_optimization (1).png" alt="topology_optimization (1)" style="zoom:80%;" />



## results

- green means linkage with tension

- red means linkage with compression

<img src="./docs/image-20210822102504200.png" alt="image-20210822102504200" style="zoom: 50%;" />

<img src="./docs/image-20210822102848149.png" alt="image-20210822102848149" style="zoom:50%;" />

<img src="./docs/image-20210822103551728.png" alt="image-20210822103551728" style="zoom:50%;" />

<img src="./docs/image-20210822103652578.png" alt="image-20210822103652578" style="zoom:50%;" />

<img src="./docs/image-20210822103841508.png" alt="image-20210822103841508" style="zoom:50%;" />

