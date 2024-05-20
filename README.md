# rayTrace

Indoor radio based ray-tracing with human blockage. 

A demo visualization.

![](README.assets/room1_pdp.gif)

Current version is only with GIF support. However, you can retrieve multipath component via `ray_matrix` and `ray_connections`. **I will reformat the code and update the functions/documents in my spare time.**

And we know that channel model for wireless sensing can be divided into two parts: target unrelated (blue and yellow rays in above GIF), and target related (red rays in above GIF) . Here, target-unrelated components are modelled via image-based method (but it could be modeled as any existing wireless communication channel)， while target-related components can be obtained from primitive-based method.

The target-related channel is still not well studied, but I have developed one novel computer-vision based approach. You can refer to my [paper](https://ieeexplore.ieee.org/document/10525191/) and another repository [testSpectrogram](https://github.com/rzy0901/testSpectrogram) for reference.

```latex
//Early access article
@ARTICLE{ren2024caster,
  author={Ren, Zhenyu and Li, Guoliang and Ji, Chenqing and Yu, Chao and Wang, Shuai and Wang, Rui},
  journal={IEEE Open Journal of the Communications Society}, 
  title={CASTER: A Computer-Vision-Assisted Wireless Channel Simulator for Gesture Recognition}, 
  year={2024},
  volume={},
  number={},
  pages={1-1},
  keywords={Videos;Wireless communication;Wireless sensor networks;Gesture recognition;Channel impulse response;Transmitters;Training;Wireless hand gesture recognition;channel model;simulation-to-reality inference},
  doi={10.1109/OJCOMS.2024.3398016}}
```

Acknowledgement: The first version of this repository was implemented on Mr. Ash Bellett's code three years ago at 2020 (When I was an undergraduate student at SUSTech). Big thanks for him! However, I could not find this codes' reference now, and if any one could find it, I will add it into this repository's citation.