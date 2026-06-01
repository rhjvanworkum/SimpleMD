"""Pygame trajectory viewer.

Replays a saved trajectory JSON file (as produced by the CLI's ``--output``
flag or the bundled ``output/`` samples). Requires the optional ``viz`` extra::

    uv sync --extra viz
    uv run python -m simplemd.viz output/test_2d_10_cd.json

The last entry of each frame is interpreted as the simulation region, following
the original viewer's convention.
"""

from __future__ import annotations

import argparse
import json
import sys

import numpy as np

BLACK = (0, 0, 0)
WHITE = (255, 255, 255)
BLUE = (0, 0, 255)
SCREEN_SIZE = 700
FPS = 60


def load_trajectory(path: str) -> np.ndarray:
    with open(path) as f:
        return np.array(json.load(f))


def run_viewer(path: str) -> None:
    """Open a window and replay the trajectory at ``path``."""
    import pygame  # imported lazily so the package imports without the viz extra

    traj = load_trajectory(path)
    pygame.init()

    loop = {"playing": False, "index": 0}
    region = traj[0][-1]
    n_atoms = len(traj[0]) - 1
    scale = SCREEN_SIZE / region[0]

    screen = pygame.display.set_mode((SCREEN_SIZE, SCREEN_SIZE), 0, 32)
    pygame.display.set_caption("MD trajectory")
    clock = pygame.time.Clock()
    font = pygame.font.Font(pygame.font.get_default_font(), 20)

    def button(text, pos, width, height, action):
        mouse = pygame.mouse.get_pos()
        click = pygame.mouse.get_pressed()
        hovering = pos[0] + width > mouse[0] > pos[0] and pos[1] + height > mouse[1] > pos[1]
        if hovering and click[0] == 1 and action is not None:
            action()
        pygame.draw.rect(screen, WHITE, (pos[0], pos[1], width, height))
        surf = font.render(text, True, BLACK)
        rect = surf.get_rect(center=(pos[0] + width / 2, pos[1] + height / 2))
        screen.blit(surf, rect)

    while True:
        screen.fill(WHITE)
        w, h = 50, 25
        button("start", [SCREEN_SIZE - w, 0], w, h, lambda: loop.update(playing=True))
        button("stop", [SCREEN_SIZE - 2 * w, 0], w, h, lambda: loop.update(playing=False))
        button("reset", [SCREEN_SIZE - 3 * w, 0], w, h, lambda: loop.update(playing=False, index=0))

        if loop["playing"] and loop["index"] < len(traj) - 1:
            loop["index"] += 1

        for i in range(n_atoms):
            pos = scale * (traj[loop["index"]][i][:2] + 0.5 * region[:2])
            pygame.draw.circle(screen, BLUE, pos, 5, 0)

        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                sys.exit()

        pygame.display.update()
        clock.tick(FPS)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Replay a SimpleMD trajectory JSON file.")
    parser.add_argument("trajectory", help="path to a trajectory JSON file")
    args = parser.parse_args(argv)
    run_viewer(args.trajectory)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
