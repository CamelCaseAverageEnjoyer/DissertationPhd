"""Проверка на работоспособность установленных библиотек. Среди них:
- kiamastro
    Нет FKIAMToolbox, не может ничего запустить
- nasapy
    Есть класс Nasa, что-то делает
- celmech
    Большая библиотека. Обязательно изучить, когда займёшься уравнениями
- poliastro
    Атмосфера: https://docs.poliastro.space/en/stable/examples/Atmospheric%20models.html
- astropy
"""

import numpy as np
import plotly.graph_objects as go

spherical_earth_map = np.load('data/map_sphere.npy')
xm, ym, zm = spherical_earth_map.T

def degree2radians(degree):
    return degree * np.pi / 180

def mapping_map_to_sphere(lon, lat, radius=1):
    # lon in [-pi, pi), lat in [-pi/2,pi/2]
    # this function maps the points of coords (lon, lat) to points onto the  sphere of radius=1

    lon = np.array(lon, dtype=np.float64)
    lat = np.array(lat, dtype=np.float64)
    lon = degree2radians(lon)
    lat = degree2radians(lat)
    xs = radius * np.cos(lon) * np.cos(lat)
    ys = radius * np.sin(lon) * np.cos(lat)
    zs = radius * np.sin(lat)
    return xs, ys, zs

fig = go.Figure(go.Scatter3d(x=xm, y=ym, z=zm, mode='lines'))
fig.update_layout(width=900, height=900)
fig.show()

