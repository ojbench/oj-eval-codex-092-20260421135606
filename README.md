# Problem 092 - 黑洞余弦的电阻网络

**ACMOJ Problem ID**: 2238

## Table of Contents

- [Problem 092 - 黑洞余弦的电阻网络](#problem-092-黑洞余弦的电阻网络)
  - [Table of Contents](#table-of-contents)
  - [Introduction](#introduction)
    - [Background](#background)
  - [Assignment Description](#assignment-description)
    - [Grade Composition](#grade-composition)
  - [Assignment Requirements](#assignment-requirements)
    - [Input Format](#input-format)
    - [Output Format](#output-format)
    - [Samples](#samples)
    - [Data Constraints](#data-constraints)
  - [Per-Testcase Resource Limits](#per-testcase-resource-limits)
  - [Test Data](#test-data)
  - [Submission Requirements](#submission-requirements)
    - [Submission Format](#submission-format)
    - [Submission Guidelines](#submission-guidelines)
    - [Evaluation Notes](#evaluation-notes)
    - [Academic Integrity](#academic-integrity)

## Introduction

### Background

This problem tests algorithmic problem-solving skills.

## Assignment Description

### Grade Composition

| Grading Component | Percentage |
| :--: | :--: |
| Pass **2238. 黑洞余弦的电阻网络** | 80% |
| Code Review | 20% |

Here are several points that need clarification:

- In the Code Review, we will **strictly examine your code style and repository organization structure, etc.**

- This assignment provides some sample data for testing, stored in the `/workspace/data/092/data_test/` directory. Note that these are not the test cases on the Online Judge. Passing all local test cases does not guarantee that you will pass the OJ tests.

- Besides the provided sample data, we also encourage you to design your own test data based on your program logic to assist debugging.

## Assignment Requirements

### Problem Description

#### 问题描述

请在 [attachment]2238材料.zip[/attachment] / [jbox](https://jbox.sjtu.edu.cn/l/S1V62p) 下载本题的下发材料。

黑洞余弦被一道有关电阻网络的物理题目难倒了。悲愤的他打算用万能的计算机来报复出题人。


设在电阻网络中：

* 有着 $n$ 个节点,编号为 $1\sim n$ ;
* 在不同节点间,有着 $m$ 条接线,编号为 $1\sim m$ ;
* 设第$k$条接线连接节点 $from(k)$ 和节点 $to(k)$ ,且 $from(k)<to(k)$ ,两个节点之间至多存在一条接线 ;
* 设节点 $k$ 上的电压为 $u_k$ ,记 $U = [u_1,...,u_n]^T$ ;
* 设节点 $k$ 上的注入电流为 $i_k$ (流出为负) ,记 $I = [i_1,...,i_n]^T$ ;
* 设接线 $k$ 上的电阻为 $r_k$ 。



黑洞余弦希望对电阻网络计算以下的特征量：

* 两节点间的等效电阻 ;
* 在给定节点电流 $I$、$u_n=0$ 的条件下，某节点的电压 ;
* 在给定节点电压 $U$ 的条件下，电阻网络的功率 。


#### 附加公式
为了帮助你计算，黑洞余弦向你提供了一些计算公式。具体见**下发文件 reading_materials.pdf**。


#### 题目要求
现在请基于以上材料，实现一个电阻网络计算器，计算材料中提到的三种特征量。善良的黑洞余弦已经帮你实现了分数类 fraction，你可以调用它的任意接口，但是请勿修改。同时他把代码模板 src.hpp 送给了你，具体均见**下发文件 fraction.hpp, src.hpp**。

你需要做的是按照要求实现电阻网络类 resistive_network 的构造函数和三个特征量计算函数。

【tips1: 测评中仅会调用 resistive_network 类的接口。对于每个测试点，resistive_network 类的构造函数仅会被调用一次，但是其他接口可能会被多次调用。】

【tips2: 为了帮助你实现 resistive_network 类的各个接口，黑洞余弦向你提供了矩阵类 matrix 的框架 (在 src.hpp 中)，你可以填补完 matrix 类并借助它来帮你实现后续内容。当然，你也可以抛开 matrix 类自己实现 resistive_network 类的各个接口，此时请将 matrix 类注释掉。】

【tips3: **除了标明可以修改的部分，请勿修改任何其他模板内容**。】

### Input Format

请按上述要求完成下发文件中的 src.hpp 并提交。

### Output Format



### Samples

No sample data provided for this problem.

### Data Constraints

| 测试点  |  描述 | 分数  |
| ------------ | ------------ | ------------ |
|  1-3 |  仅调用 `get_power` 函数，不检查内存泄漏 |  21 |
|  4 | 仅调用 `get_equivalent_resistance` 函数，$n,m\le 5$，不检查内存泄漏  | 7  |
|  5-6 |  仅调用 `get_equivalent_resistance` 函数，$n\le 10$，不检查内存泄漏 | 14  |
|  7-8 |  仅调用 `get_voltage` 函数，$n\le 10$，不检查内存泄漏 | 14  |
|  9-10 |  resistive_network 类的所有函数均会调用，$n\le 10$，不检查内存泄漏 | 14  |
| 11-20 | 检查内存泄漏，其余与1-10分别相同 | 30 |

## Per-Testcase Resource Limits

- **Time Limit (per test case)**: 1500 ms
- **Memory Limit (per test case)**: 512 MiB
- **Disk Usage**: No disk usage is permitted.

## Test Data

The test data for this problem is located at `/workspace/data/092/data_test/`.

Attachments and starter files for this problem are located at `/workspace/data/092/`.

## Submission Requirements

### Submission Format

Submit a single C++ source file to the Online Judge. The OJ provides its own `main.cpp` which `#include "src.hpp"`. Your submission will be written to `src.hpp` and compiled together with the OJ's driver code.

### Submission Guidelines

- The submitted code must be able to compile successfully through the above compilation process
- The program needs to be able to read data from standard input and write results to standard output
- Please ensure the code runs correctly within the given time and space limits
- **You must use C++ or C language** to implement this assignment

### Evaluation Notes

- The evaluation system will test your program using the provided test data
- The program output must exactly match the expected output (including format)
- Exceeding time or memory limits will be judged as the corresponding error type
- Please pay attention to the overall time performance of your code and the time complexity of each part of your algorithm.

### Academic Integrity

If any violations are found during evaluation or code review (including but not limited to using unconventional methods to pass test cases), your final score may be significantly reduced or become **0 points**.
