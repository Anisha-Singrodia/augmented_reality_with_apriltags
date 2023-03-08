# CIS580 HW2 

To run this program:

```
cd code
python main.py
```

We also provided some helper flags. Please check `main.py` for details. You can generate your visualizations with either PnP or P3P algorithm, but you still need to implement both of them. Although you are not asked to implement the renderer, you are still encouraged to look through the code of renderer as you may need to renderer your results by yourself next time. 

PS: remember to complete the `est_homography.py` with the function you just wrote for HW1.

## Debugging

It's recommended to run the program with `--debug` when you start to work on this homework since the rendering takes about 2 mins to finish on a PC. 

Also, note that the main program has several other args you can set, please have a look at line 40 in the `main.py` for more details: You can pass `--solver PnP` or `--solver P3P` to the program to toggle between the solving methods.

Note, we also provide the `.vscode` launch configuration for you to easily debug in VSCode.

## Customization

You need to assign different values to  `click_point` in `main.py` to render the drill at different places. 



# augmented_reality_with_apriltags
