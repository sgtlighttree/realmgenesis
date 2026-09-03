export const shouldUseGLFill = (o: {
  styleId: string; projectionType: string; webglAvailable: boolean; contextLost: boolean;
}): boolean =>
  o.styleId === 'default' &&
  o.projectionType !== 'dymaxion' &&
  o.webglAvailable &&
  !o.contextLost;
