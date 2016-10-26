import pygame
import random
import math
import numpy as np
from typing import List

background_colour = (0,0,0)
(width, height) = (400, 400)

class Particle():
    def __init__(self):
        self.size = 20
        self.sizesq = self.size*self.size
        self.x = random.randint(self.size, width-self.size)
        self.y = random.random()*width
        self.colour = (random.randint(0,255), random.randint(0,255), random.randint(0,255))
        self.thickness = 0 #The width is set to 0 in order to get a filled shape
        self.speed = 15
        self.angle = random.uniform(0, math.pi*2)

    def display(self):
        pygame.draw.circle(screen, self.colour, (int(self.x), int(self.y)), self.size, self.thickness)

    def move(self):
        self.x += math.sin(self.angle) * self.speed
        self.y -= math.cos(self.angle) * self.speed

    def bounce(self):
        if self.x > width - self.size:
            self.x = 2*(width - self.size) - self.x
            self.angle = - self.angle

        elif self.x < self.size:
            self.x = 2*self.size - self.x
            self.angle = - self.angle

        if self.y > height - self.size:
            self.y = 2*(height - self.size) - self.y
            self.angle = math.pi - self.angle

        elif self.y < self.size:
            self.y = 2*self.size - self.y
            self.angle = math.pi - self.angle

def CheckColide(a, b):
    """
    :param a: Particle
    :param b: Particle
    :return: bool
    """
    return (b.x-a.x)*(b.x-a.x)+(b.y-a.y)*(b.y-a.y)<(a.sizesq+b.sizesq)

def Collide(partlist):
    """

    :param partlist: List[Particle]
    :return:
    """
    npars=len(partlist)

    #Build Collision Matrix, so far it does not include multiple collisions
    coll_mat=np.zeros((npars,npars),dtype=bool)
    rows = range(1,npars)
    for i in rows:
        for j in range(i):
            coll_mat[i,j]=CheckColide(partlist[i],partlist[j])
    for i in rows:
        for j in range(i):
            if coll_mat[i,j]:
                old_ang=partlist[i].angle
                partlist[i].angle=partlist[j].angle
                partlist[j].angle =old_ang
                print "hi"
    # for i in partlist:
    #     for j in partlist:
    #         if CheckColide(i,j):
    #             old_ang = i.angle
    #             i.angle=j.angle
    #             i.angle =old_ang
    #             print "hi"




screen = pygame.display.set_mode((width, height))
pygame.display.set_caption('Bouncin')

number_of_particles = 15
my_particles = []

for n in range(number_of_particles):
    my_particles.append(Particle())

running=True

while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    screen.fill(background_colour)
    Collide(my_particles)
    for particle in my_particles:
        particle.move()
        particle.bounce()
        particle.display()

    pygame.display.flip()