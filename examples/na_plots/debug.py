xhat = (0.16218229589602626, 0.59900602672105474) # 1760
x0 = (0.1434482564694205, 0.57756856005299395) # 1761
x1 = (0.14691397003969292, 0.5493104985513505)
x2 = (0.11498063497341973, 0.57721679771522882) # 1759
x3 = (0.12921249311229097, 0.60832241498586515)
x4 = (0.20403924632836226, 0.56817158177670701) # 1748
x5 = (0.19025455987981871, 0.59328631016751254) # 1777

DT0 = (0.24072380496331644, 0.97059365839880862)
DT1 = (0.25802772770183324, 0.96613751181549135)
DT2 = (0.19507410952437368, 0.98078850512904802)
DT3 = (0.20749124394337137, 0.97823687503939049)
DT4 = (0.33769745963727349, 0.94125470822436375)
DT5 = (0.30502789008552139, 0.95234341824258706)

points = [(0, 0), x0, x1, xhat]

xmin = min(_[0] for _ in points)
xmax = max(_[0] for _ in points)
ymin = min(_[1] for _ in points)
ymax = max(_[1] for _ in points)

w, h = xmax - xmin, ymax - ymin
xmin -= 0.1*max(w, h)
xmax += 0.1*max(w, h)
ymin -= 0.1*max(w, h)
ymax += 0.1*max(w, h)

s, scale, width = 15, 10, 0.01

plt.figure()
plt.plot([0, xhat[0]], [0, xhat[1]], linestyle='--', c='k', linewidth=1, zorder=1)
for z in [xhat, x0, x1]:
    plt.plot([0, z[0]], [0, z[1]], linestyle='-', c='gray', linewidth=0.5, zorder=1)
plt.scatter(0, 0, marker='x', c='k')
plt.scatter(*xhat, facecolor='cyan', edgecolor='k', s=s, zorder=2)
plt.scatter(*x0, facecolor='red', edgecolor='k', s=s, zorder=2)
plt.scatter(*x1, facecolor='orange', edgecolor='k', s=s, zorder=2)
plt.scatter(*x2, facecolor='pink', edgecolor='k', s=s, zorder=2)
plt.scatter(*x3, facecolor='brown', edgecolor='k', s=s, zorder=2)
plt.scatter(*x4, facecolor='yellow', edgecolor='k', s=s, zorder=2)
plt.scatter(*x5, facecolor='purple', edgecolor='k', s=s, zorder=2)
# plt.scatter(*x6, facecolor='blue', edgecolor='k', s=s, zorder=2)
# plt.scatter(*xlam, c='pink')
plt.quiver(*x0, *DT0, color='black', scale=scale, width=width, zorder=1)
plt.quiver(*x1, *DT1, color='black', scale=scale, width=width, zorder=1)
plt.quiver(*x2, *DT2, color='black', scale=scale, width=width, zorder=1)
plt.quiver(*x3, *DT3, color='black', scale=scale, width=width, zorder=1)
plt.quiver(*x4, *DT4, color='black', scale=scale, width=width, zorder=1)
plt.quiver(*x5, *DT5, color='black', scale=scale, width=width, zorder=1)
# plt.quiver(*xhat, *DT, color='cyan')
plt.plot([-1, 1, 1, -1, -1], [-1, -1, 1, 1, -1], c='k', linewidth=2)
plt.gca().set_aspect('equal')
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.tight_layout()
