import turtle
import os


class domainpositonkind:
    r'''domainkind 0:cs,1:ds,2:TMhelix'''
    def __init__(self, positionstart, positionend, domainkind):
        self.start = positionstart
        self.end = positionend
        self.kind =domainkind

def drawdomain(domanlist):
    turtle.setup(600, 30, None, None)
    turtle.pensize(1)
    turtle.hideturtle()
    height = 1
    width = 10
    turtle.addshape("bar", ((width / 2, 0), (-width / 2, 0), (-width / 2, height), (width / 2, height)))
    turtle.shape('bar')
   # turtle.pendown()
    turtle.write('this is text')
    turtle.penup()
    turtle.goto(150, width/2)
    for dk in domanlist:
        if dk.kind == 1:
            turtle.pencolor('yellow')
        if dk.kind == 2:
            turtle.pencolor('blue')
        for i in range(1, dk.end-dk.start):
            turtle.stamp()
            turtle.fd(1)
    turtle.write('this is end text')
    turtle.penup()
    ts = turtle.getscreen()
    path = os.path.abspath(os.curdir) + '\\'
    ts.getcanvas().postscript(file= path+ "duck.eps")
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