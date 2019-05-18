import turtle
import os
from PIL import Image


class domainpositonkind:
    r'''domainkind 0:cs,1:ds,2:TMhelix'''
    def __init__(self, positionstart, positionend, domainkind):
        self.start = positionstart
        self.end = positionend
        self.kind =domainkind
def saveImg():
    global canvas
    path = os.path.abspath(os.curdir) + '\\'
    canvas.postscript(file=path+"duck.eps", colormode='color')
    imgNew = Image.open(path+"duck.eps")
    imgNew.convert("RGBA")
    imgNew.thumbnail((2000, 2000), Image.ANTIALIAS)
    imgNew.save(path+'testImg.png', quality=90)


def drawdomain(domanlist):
    scrwidth = 400
    #todo 计算画布宽度
    turtle.setup(scrwidth, 140, None, None)
    turtle.pensize(1)
    turtle.hideturtle()
    height = 1
    width = 40
    turtle.addshape("bar", ((width / 2, 0), (-width / 2, 0), (-width / 2, height), (width / 2, height)))
    turtle.shape('bar')
    turtle.penup()
    turtle.goto(-scrwidth/2+40, -width/2)
    turtle.write('this is text',True, align='left', font=("Arial", 24, "normal"))
    turtle.goto(turtle.xcor(), 0)
    turtle.pendown()
    for dk in domanlist:
        if dk.kind == 1:
            turtle.pencolor('yellow')
        if dk.kind == 2:
            turtle.pencolor('blue')
        for i in range(1, dk.end-dk.start):
            turtle.stamp()
            turtle.fd(1)
    turtle.penup()
    turtle.goto(turtle.xcor(), -width/2)
    turtle.write('this is text', True, align='left', font=("Arial", 24, "normal"))
    turtle.penup()
    ts = turtle.getscreen()
    path = os.path.abspath(os.curdir) + '\\'
    ts.getcanvas().postscript(file= path+ "duck.eps", colormode='color')
    try:
        imgNew = Image.open(path + "duck.eps")
        imgNew.convert("RGBA")
        imgNew.thumbnail((1000, 200), Image.ANTIALIAS)
        imgNew.save(path + 'testImg.png', quality=70)
    except Exception as e:
        print(e)
    pass
#    turtle.done()
    turtle.bye()

if __name__ == '__main__':
    adk = domainpositonkind(1,30,1)
    cdk = domainpositonkind(45,80,2)
    ddk = domainpositonkind(90,120,1)
    edk = domainpositonkind(120,130,1)
    fdk = domainpositonkind(135,180,2)
    c = []
    c.append(adk)
    c.append(cdk)
    c.append(ddk)
    c.append(edk)
    c.append(fdk)
    drawdomain(c)
    #saveImg()

