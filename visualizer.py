import pygame,sys
from pygame.locals import *
import json
import numpy as np

# initialization
file = 'output/test_2d_10_cd.json'
with open(file) as f:
    data = json.load(f)
traj = np.array(data)

pygame.init()

# settings
screenSize = 700

loop = False

FPS = 60
frameCount = 0
updateFreq = 1
trajIndex = 0
fpsClock = pygame.time.Clock()

region = traj[0][-1]
nAtoms = len(traj[0]) - 1
scaleFactor = screenSize / region[0]

black = (  0,   0,   0)
white = (255, 255, 255)
red   = (255,   0,   0)
green = (  0, 255,   0)
blue  = (  0,   0, 255)

#set up the window
screen = pygame.display.set_mode((screenSize, screenSize), 0, 32)
pygame.display.set_caption('MD trajectory')

# add button function
def addButton(screen, text, color, pos, width, height, action):
    mouse = pygame.mouse.get_pos()
    click = pygame.mouse.get_pressed()

    if pos[0]+width > mouse[0] > pos[0] and pos[1]+height > mouse[1] > pos[1]:
        if click[0] == 1 and action != None:
            action()

    pygame.draw.rect(screen, color, (pos[0], pos[1], width, height))
    smallText = pygame.font.Font("freesansbold.ttf", 20)
    textSurf = smallText.render(text, True, black)
    textRect = textSurf.get_rect()
    textRect.center = ((pos[0] + (width / 2)), (pos[1] + (height / 2)))
    screen.blit(textSurf, textRect)

def startLoop():
    global loop
    loop = True

def stopLoop():
    global loop
    loop = False

def resetLoop():
    global trajIndex, loop
    trajIndex = 0
    loop = False

# the main game loop
while True:
    screen.fill(white)

    width = 50
    height = 25
    addButton(screen, 'start', white, [screenSize - width, 0], width, height, startLoop)
    addButton(screen, 'stop', white, [screenSize - 2*width, 0], width, height, stopLoop)
    addButton(screen, 'reset', white, [screenSize - 3*width, 0], width, height, resetLoop)

    if (loop):
        if (frameCount % updateFreq == 0):
            if(trajIndex < len(traj) - 1): trajIndex += 1

    for i in range(nAtoms):
        pygame.draw.circle(screen, blue, scaleFactor * (traj[trajIndex][i][:2] + 0.5 * region[:2]), 5, 0)

    for event in pygame.event.get():
        if event.type == QUIT:
            pygame.quit()
            sys.exit()

    frameCount += 1
    pygame.display.update()
    fpsClock.tick(FPS)