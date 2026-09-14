#!/usr/bin/env python3

"""Generate the README figures from the simulator's CSV output."""

from __future__ import annotations

import csv
import math
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "assets"
WIDTH, HEIGHT = 1600, 900
PLOT = (150, 120, 1515, 760)

INK = "#172033"
MUTED = "#596579"
GRID = "#d8deea"
BLUE = "#2378c9"
CYAN = "#10a6a6"
ORANGE = "#e9772e"
RED = "#cf3f49"


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    filename = "Arial Bold.ttf" if bold else "Arial.ttf"
    candidates = [
        Path("/System/Library/Fonts/Supplemental") / filename,
        Path("/Library/Fonts") / filename,
        Path("/usr/share/fonts/truetype/dejavu")
        / ("DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"),
    ]
    for path in candidates:
        if path.exists():
            return ImageFont.truetype(str(path), size)
    return ImageFont.load_default(size=size)


TITLE_FONT = font(34, bold=True)
BODY_FONT = font(21)
SMALL_FONT = font(18)
LABEL_FONT = font(20)


def read_csv(filename: str) -> list[dict[str, float]]:
    with (ROOT / filename).open(newline="", encoding="utf-8") as source:
        return [
            {key: float(value) for key, value in row.items()}
            for row in csv.DictReader(source)
        ]


def new_chart(title: str) -> tuple[Image.Image, ImageDraw.ImageDraw]:
    image = Image.new("RGB", (WIDTH, HEIGHT), "white")
    draw = ImageDraw.Draw(image)
    draw.text((75, 45), title, fill=INK, font=TITLE_FONT)
    return image, draw


def linear(value: float, low: float, high: float, start: float, end: float) -> float:
    return start + (value - low) / (high - low) * (end - start)


def logarithmic(value: float, low: float, high: float,
                start: float, end: float) -> float:
    return linear(math.log10(value), math.log10(low), math.log10(high), start, end)


def draw_axes(
    draw: ImageDraw.ImageDraw,
    x_ticks: list[tuple[float, str]],
    y_ticks: list[tuple[float, str]],
    x_title: str,
    y_title: str,
) -> None:
    left, top, right, bottom = PLOT
    for x, label in x_ticks:
        draw.line((x, top, x, bottom), fill=GRID, width=1)
        draw.text((x, bottom + 17), label, fill=MUTED, font=SMALL_FONT,
                  anchor="ma")
    for y, label in y_ticks:
        draw.line((left, y, right, y), fill=GRID, width=1)
        draw.text((left - 18, y), label, fill=MUTED, font=SMALL_FONT,
                  anchor="rm")
    draw.line((left, top, left, bottom), fill="#9aa5b5", width=2)
    draw.line((left, bottom, right, bottom), fill="#9aa5b5", width=2)
    draw.text(((left + right) / 2, 845), x_title, fill=MUTED,
              font=LABEL_FONT, anchor="mm")

    label = Image.new("RGBA", (650, 50), (255, 255, 255, 0))
    label_draw = ImageDraw.Draw(label)
    label_draw.text((325, 25), y_title, fill=MUTED, font=LABEL_FONT,
                    anchor="mm")
    rotated = label.rotate(90, expand=True)
    y_position = int((top + bottom - rotated.height) / 2)
    draw._image.paste(rotated, (25, y_position), rotated)


def draw_series(draw: ImageDraw.ImageDraw, points: list[tuple[float, float]],
                colour: str, width: int = 5, radius: int = 7) -> None:
    draw.line(points, fill=colour, width=width, joint="curve")
    for x, y in points:
        draw.ellipse((x - radius, y - radius, x + radius, y + radius),
                     fill=colour, outline="white", width=2)


def draw_legend(draw: ImageDraw.ImageDraw,
                entries: list[tuple[str, str]]) -> None:
    x, y = 175, 95
    for label, colour in entries:
        draw.line((x, y, x + 42, y), fill=colour, width=6)
        draw.text((x + 55, y), label, fill=INK, font=BODY_FONT, anchor="lm")
        x += 250


def plot_benchmark() -> None:
    rows = read_csv("benchmark.csv")
    x_min, x_max = 80.0, 12500.0
    y_min, y_max = 0.03, 160.0
    left, top, right, bottom = PLOT
    map_x = lambda value: logarithmic(value, x_min, x_max, left, right)
    map_y = lambda value: logarithmic(value, y_min, y_max, bottom, top)

    image, draw = new_chart("Force-calculation runtime")
    counts = [row["particle_count"] for row in rows]
    x_ticks = [(map_x(value), f"{int(value):,}") for value in counts]
    y_values = [0.05, 0.1, 0.5, 1, 5, 10, 50, 100]
    y_ticks = [(map_y(value), f"{value:g}") for value in y_values]
    draw_axes(draw, x_ticks, y_ticks, "Particle count (log scale)",
              "Median runtime (ms, log scale)")

    direct = [(map_x(row["particle_count"]), map_y(row["direct_ms"]))
              for row in rows]
    barnes_hut = [(map_x(row["particle_count"]), map_y(row["barnes_hut_ms"]))
                  for row in rows]
    draw_series(draw, direct, ORANGE)
    draw_series(draw, barnes_hut, CYAN)
    draw_legend(draw, [("Direct O(n²)", ORANGE), ("Barnes–Hut", CYAN)])

    final_x, final_y = barnes_hut[-1]
    draw.text((final_x - 20, final_y + 42),
              f"{rows[-1]['speedup']:.2f}× faster", fill=INK,
              font=BODY_FONT, anchor="ra")
    image.save(ASSETS / "benchmark-scaling.png", optimize=True)


def plot_energy_drift() -> None:
    rows = read_csv("integrator_accuracy.csv")
    x_min, x_max = 0.0, rows[-1]["time"]
    y_min, y_max = 1e-12, 10.0
    left, top, right, bottom = PLOT
    map_x = lambda value: linear(value, x_min, x_max, left, right)
    map_y = lambda value: logarithmic(max(value, y_min), y_min, y_max,
                                      bottom, top)

    image, draw = new_chart("Relative energy drift in a two-body orbit")
    x_values = [0, 5, 10, 15, 20]
    x_ticks = [(map_x(value), f"{value:g}") for value in x_values]
    y_values = [1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1, 10]
    y_ticks = [(map_y(value), f"{value:g}") for value in y_values]
    draw_axes(draw, x_ticks, y_ticks, "Simulated time",
              "Relative energy drift (%) — log scale")

    euler = [(map_x(row["time"]), map_y(row["euler_energy_drift"] * 100.0))
             for row in rows]
    verlet = [(map_x(row["time"]), map_y(row["verlet_energy_drift"] * 100.0))
              for row in rows]
    draw.line(euler, fill=RED, width=5, joint="curve")
    draw.line(verlet, fill=BLUE, width=5, joint="curve")
    draw_legend(draw, [("Euler", RED), ("Velocity Verlet", BLUE)])
    draw.text((right - 15, euler[-1][1] + 30),
              f"Euler: {rows[-1]['euler_energy_drift'] * 100:.2f}%",
              fill=RED, font=SMALL_FONT, anchor="ra")
    draw.text((right - 15, verlet[-1][1] - 28),
              f"Verlet: {rows[-1]['verlet_energy_drift'] * 100:.2e}%",
              fill=BLUE, font=SMALL_FONT, anchor="ra")
    image.save(ASSETS / "integrator-energy-drift.png", optimize=True)


def plot_theta_tradeoff() -> None:
    rows = read_csv("theta_accuracy.csv")
    errors = [row["relative_acceleration_error"] * 100.0 for row in rows]
    speedups = [row["speedup"] for row in rows]
    x_min, x_max = 0.0, 2.1
    y_min, y_max = 0.0, 18.0
    left, top, right, bottom = PLOT
    map_x = lambda value: linear(value, x_min, x_max, left, right)
    map_y = lambda value: linear(value, y_min, y_max, bottom, top)

    image, draw = new_chart("Barnes–Hut speed–accuracy trade-off")
    x_values = [0, 0.5, 1.0, 1.5, 2.0]
    y_values = [0, 4, 8, 12, 16]
    x_ticks = [(map_x(value), f"{value:g}") for value in x_values]
    y_ticks = [(map_y(value), f"{value:g}×") for value in y_values]
    draw_axes(draw, x_ticks, y_ticks, "Relative acceleration error (%)",
              "Speedup over direct calculation")

    points = [(map_x(error), map_y(speedup))
              for error, speedup in zip(errors, speedups)]
    draw_series(draw, points, BLUE, width=4, radius=9)
    for (x, y), row in zip(points, rows):
        draw.text((x + 13, y - 13), f"θ = {row['theta']:g}", fill=INK,
                  font=SMALL_FONT, anchor="ls")
    image.save(ASSETS / "theta-tradeoff.png", optimize=True)


def main() -> None:
    ASSETS.mkdir(parents=True, exist_ok=True)
    plot_benchmark()
    plot_energy_drift()
    plot_theta_tradeoff()
    print("Generated README figures in assets/")


if __name__ == "__main__":
    main()
