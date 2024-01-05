# Passive Handwriting Tracking

This repository is the course Project in the course CSE5010: Wireless Network and Mobile Computing in 2023 Fall by lecturer Prof. Jin Zhang. 

![HitCount](https://img.shields.io/endpoint?url=https%3A%2F%2Fhits.dwyl.com%2FJcq242818%2FPassive-Handwriting-Tracking.json%3Fcolor%3Dpink)
![Static Badge](https://img.shields.io/badge/Matlab-2023a-salmon)
![Static Badge](https://img.shields.io/badge/Python-3.9.13-3776AB?logo=python)
![Static Badge](https://img.shields.io/badge/OpenCV-4.7.0.72-5C3EE8?logo=opencv)
![Static Badge](https://img.shields.io/badge/Ubuntu-18.04-red?logo=ubuntu)

## How to Use

### 下载程序
   
下载`scripts/scripts_jcq`文件夹下所有文件以及`scripts/scripts_hyx/frameAnalyse.ipynb`到本地。


### 下载数据

下载`data`文件夹下所有文件到本地。

**Note:** The dataset in the `./Data/Broken_data` folder is the wrong data containing noise collected due to the **device clock problem** when we test the tracking results under the condition that the reference channel was NLOS in the last week. In order to prove the authenticity of the Figure on the last page of our PPT, this data is specially retained. The TAs can also change the path to this folder to run the results.


### 修改程序

每次运行时，需修改如下几个部分：

#### main_jcq.m

**Firstly**, Please open the project folder，the `main.m` file contain all the algorithm included in this project. So, please run the `main.m` file to get the results mentioned in the PPT directly.

Note that in the `main.m` file, there are examples as follows:

```matlab
filename_data = "E:\Desktop\Project\Data\Right_data\3.bin";
```

So if you want to run predictions for different handwritten number tracks, you need to **modify the file path for this line of in the main.m file**. All the collected handwritten numerals dataset are placed in the `./Data/Right_data` folder. All you need to do to run this file is to change the corresponding dataset path.


#### frameAnalyse.ipynb

2块2行，文件路径：

```python
dataPath = r'D:\Github\Passive-Handwriting-Tracking\data'
```

2块3行，研究对象：`zero`、`two`、`three`、`seven`、`eight`、`star`。

```python
shapeType = 'star'
```

### 运行程序

#### main_jcq.m

一键运行 -> 弹窗查看结果

#### frameAnalyse.ipynb

重启内核 -> 一键运行所有 -> 底部查看结果

**注意：**请尽量确保以下版本匹配：

- Python：3.9.13
- OpenCV：4.7.0.72
- SciPy：1.9.1
- NumPy：1.21.5
- statsmodels：0.13.2


## Citation

The paper based on this project (see in the below citation) has been submitted to IEEE WCL. If you find the repository is helpful to your project, please cite as follows:

```latex
@misc{yu2023passive,
      title={Passive Handwriting Tracking via Weak mmWave Communication Signals},
      author={Chao Yu and Yan Luo and Renqi Chen and Rui Wang},
      year={2023},
      eprint={2311.01781},
      archivePrefix={arXiv},
      primaryClass={eess.SY}
}
```

## Reference

- Kun Qian1, Chenshu Wu2, Yi Zhang1, Guidong Zhang3, Zheng Yang1, Yunhao Liu1,4. 2018. Widar2.0: Passive Human Tracking with a Single Wi-Fi Link. In MobiSys ’18: The 16th Annual International Conference on Mobile Systems, Applications, and Services, June 10–15, 2018, Munich, Germany. ACM, New York, NY, USA, 12 pages. https://doi.org/10.1145/3210240.3210314


