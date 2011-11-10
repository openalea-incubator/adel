color_list=[(0,0,0),
            (255,0,0),
            (0,255,0),
            (0,0,255),
            (255,255,0),
            (0,255,255),
            (255,0,255),
            (128,255,0),
            (0,128,255),
            (255,0,128),
            (0,255,128),
            (128,0,255),
            (255,128,0),
            (128,128,255),
            (255,128,128),
            (128,255,128),
            (255,255,255)
            ]

def col_item (ind, color_list=color_list) :
    if ind is None:
        return lambda x: color_list[(x-1) % len(color_list)]
    else:
        return color_list[(ind-1) % len(color_list)],

