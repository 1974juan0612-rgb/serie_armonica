# Despliegue en GitHub Pages

## Requisito
Una cuenta de GitHub (gratuita). Si no tienes: https://github.com/signup

## Paso 1 — Crear el repositorio en GitHub

1. Ve a https://github.com/new
2. Nombre del repositorio: `serie_armonica` (o el que quieras)
3. Déjalo en **público**
4. **NO** marques "Add a README", "Add .gitignore" ni "Choose a license"
5. Haz clic en **Create repository**

## Paso 2 — Subir el código

En la página siguiente, GitHub te mostrará comandos. Copia y pega estos en tu terminal:

```bash
git remote add origin https://github.com/TU-USUARIO/serie_armonica.git
git branch -M master
git push -u origin master
```

Reemplaza `TU-USUARIO` por tu nombre de usuario de GitHub.

## Paso 3 — Activar GitHub Pages

1. Ve a https://github.com/TU-USUARIO/serie_armonica/settings/pages
2. En **Source**, selecciona **Deploy from a branch**
3. Branch: `master`, carpeta: `/ (root)`
4. Guarda

En 1-2 minutos tu sitio estará vivo en:
```
https://TU-USUARIO.github.io/serie_armonica/
```

## Sitio publicado

| Página | URL |
|--------|-----|
| Portal principal | `https://TU-USUARIO.github.io/serie_armonica/` |
| Dashboard CERN | `https://TU-USUARIO.github.io/serie_armonica/dashboard.html` |
| Malla Elástica | `https://TU-USUARIO.github.io/serie_armonica/interfaz/malla_musculo.html` |

## Notas

- El dashboard conecta en vivo con la API de HepData (CERN). Si falla la conexión,
  usa automáticamente una rampa de presión estandarizada local.
- Los experimentos Python requieren ejecución local (no corren en el navegador).
- Para añadir un dominio personalizado, configúralo en Settings > Pages > Custom domain.
