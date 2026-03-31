library(ggplot2)
library(hexSticker)
library(showtext)
library(ggtext)

font_add_google("Fira Sans", "fira")
showtext_auto()

# --- Build the avocado-as-censored-distribution plot ---

# Shift everything right by `pad` so the left edge of the avocado (the
# censoring wall) isn't at the very edge of the plotting region.
pad <- 1.5

# Avocado shape: right-skewed gamma density, censored at x = pad
x_raw <- seq(0, 5, length.out = 300)
y_raw <- dgamma(x_raw, shape = 2.2, rate = 1.1)
y_raw <- y_raw / max(y_raw) * 3.4
x_shifted <- x_raw + pad

# Outer avocado polygon (dark green skin)
avocado_df <- data.frame(
  x = c(pad, x_shifted, max(x_shifted), pad),
  y = c(0, y_raw, 0, 0)
)

# Inner flesh — lighter, inset from skin
inset <- 0.14
x_inner <- seq(inset, 4.8, length.out = 300) + pad
y_inner <- dgamma(seq(inset, 4.8, length.out = 300), shape = 2.2, rate = 1.1)
y_inner <- y_inner / max(dgamma(x_raw, shape = 2.2, rate = 1.1)) * 3.4 * 0.87
inner_df <- data.frame(
  x = c(pad + inset, x_inner, max(x_inner), pad + inset),
  y = c(inset, y_inner, inset, inset)
)

# Skin band on the outer curve
skin_y_outer <- dgamma(x_raw, shape = 2.2, rate = 1.1)
skin_y_outer <- skin_y_outer / max(skin_y_outer) * 3.4
skin_y_inner <- skin_y_outer * 0.90
skin_df <- data.frame(
  x = c(x_shifted, rev(x_shifted)),
  y = c(skin_y_outer, rev(skin_y_inner))
)

# --- The pit: full circle at the centroid of the distribution ---
pit_r <- 0.65
theta <- seq(0, 2 * pi, length.out = 200)

# Centroid of Gamma(2.2, 1.1): mean = shape/rate = 2.0
gamma_shape <- 2.2
gamma_rate <- 1.1
x_mean <- gamma_shape / gamma_rate  # = 2.0

# Compute the actual 2D centroid of the filled area under the curve
# y_bar = integral(f(x)^2 dx) / (2 * integral(f(x) dx))
# For visual placement: use ~40% of peak height at the x-centroid
peak_height <- max(y_raw)  # = 3.4 (already scaled)
pit_cx <- pad + x_mean - 0.45
pit_cy <- peak_height * 0.35
pit_df <- data.frame(
  x = pit_r * cos(theta) + pit_cx,
  y = pit_r * sin(theta) + pit_cy
)

# Pit shadow ring
sr <- pit_r + 0.05
shadow_df <- data.frame(
  x = sr * cos(theta) + pit_cx,
  y = sr * sin(theta) + pit_cy
)

# Pit highlight (3D sheen)
hr <- 0.20
theta_h <- seq(0, 2 * pi, length.out = 100)
highlight_df <- data.frame(
  x = hr * cos(theta_h) + pit_cx + 0.10,
  y = hr * sin(theta_h) + pit_cy + 0.08
)

# --- Assemble the subplot ---
p <- ggplot() +
  # Dark green skin (full shape)
  geom_polygon(data = avocado_df, aes(x, y),
               fill = "#4A6E1E", color = NA) +
  # Lighter flesh interior
  geom_polygon(data = inner_df, aes(x, y),
               fill = "#C5D99B", color = NA) +
  # Darker skin band on curve
  geom_polygon(data = skin_df, aes(x, y),
               fill = "#3D5C18", color = NA) +
  # Left edge skin strip
  annotate("rect", xmin = pad - 0.02, xmax = pad + inset,
           ymin = 0, ymax = max(y_raw) * 0.99, fill = "#4A6E1E") +
  # Bottom skin strip
  annotate("rect", xmin = pad, xmax = pad + 4.5,
           ymin = 0, ymax = inset, fill = "#4A6E1E") +
  # Pit shadow
  geom_polygon(data = shadow_df, aes(x, y),
               fill = "#3A2510", color = NA) +
  # Pit
  geom_polygon(data = pit_df, aes(x, y),
               fill = "#6B3A1F", color = "#4A2810", linewidth = 0.4) +
  # Pit highlight
  geom_polygon(data = highlight_df, aes(x, y),
               fill = "#8B5E3C", color = NA, alpha = 0.7) +
  # Y-axis at x = pad (the censoring wall)
  geom_segment(aes(x = pad, xend = pad, y = 0, yend = 3.5),
               color = "white", linewidth = 0.45) +
  # X-axis
  geom_segment(aes(x = pad, xend = pad + 5.2, y = 0, yend = 0),
               color = "white", linewidth = 0.45) +
  coord_cartesian(xlim = c(0, 7.3), ylim = c(-0.15, 3.9)) +
  theme_void() +
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA))

# --- Create hex sticker ---
s <- sticker(
  p,
  package = "",
  s_x = 1.10,
  s_y = 1.15,
  s_width = 1.85,
  s_height = 1.45,
  h_fill = "#2C3E50",
  h_color = "#ECF0F1",
  h_size = 1.6,
  url = "",
  filename = "/Users/hsantanna/repos/bunching/inst/logo/logo.png",
  dpi = 300
)

# Add "bRunching" with green R
s <- s +
  geom_richtext(
    aes(x = 1, y = 0.40),
    label = "<span style='color:#FFFFFF;'>b</span><span style='color:#7AB648;'>R</span><span style='color:#FFFFFF;'>unching</span>",
    size = 7.5, family = "fira", fontface = "bold",
    fill = NA, label.color = NA
  )

ggsave("/Users/hsantanna/repos/bunching/inst/logo/logo.png", s,
       width = 43.9, height = 50.8, units = "mm", dpi = 300,
       bg = "transparent")

message("Logo saved to inst/logo/logo.png")
