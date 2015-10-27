# from terminal in this directory, run py.test

from bumps import crop
from numpy import ones, sum, mgrid
from numpy.random import rand

# first simple test image is a rect of 0's on a background of 1's
rect = ones([200,200])
rect = rect>0
rect[50:120,50:130] = False

# second simple test image is a circle of 1's on a bg of 0's
rows, cols = mgrid[-100:100,-100:100]
circle = (rows**2+cols**2) < 650

def test_crop_size_on_rect():
    res  = crop(rect)
    assert sum(res[0])==0

def test_crop_size_on_circle():
    res  = crop(circle)
    assert sum(res[0])==2029

def test_crop_size_on_rect_with_padding():
    res  = crop(rect,4)
    assert sum(res[0])==1264

def test_crop_location_on_rect():
    res  = crop(rect)
    assert res[1]==[50,120,50,130]

def test_crop_location_on_rectangle():
    random_rect = ones([300,200])
    random_rect = random_rect>0
    a = int(rand()*300)
    b = int(rand()*(300-(a+1)))+a+1

    c = int(rand()*200)
    d = int(rand()*(200-(c+1)))+c+1

    print [a,b,c,d]

    random_rect[a:b,c:d] = False
    res  = crop(random_rect)
    assert res[1]==[a,b,c,d]