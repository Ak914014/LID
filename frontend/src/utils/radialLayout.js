/**
 * Professional Radial Layout Algorithm for Protein-Ligand Interaction Diagrams
 * Similar to LigPlot/Maestro style layout
 * 
 * Features:
 * - Deterministic geometry-based layout
 * - Preserves biological orientation (3D angles)
 * - Collision detection and resolution
 * - Curved Bézier connectors
 * - Residue type grouping
 */

/**
 * Calculate bounding box for a residue label
 */
export function getResidueBoundingBox(residue, radius) {
  const labelWidth = 50; // Approximate width for "XXX A: 999"
  const labelHeight = 30; // Height for two-line label
  const padding = 8;
  
  return {
    width: labelWidth + padding * 2,
    height: labelHeight + padding * 2,
    radius: radius,
    angle: residue.angle || 0,
  };
}

/**
 * Check if two bounding boxes overlap
 */
export function boxesOverlap(box1, box2, minSeparation = 5) {
  // Convert to cartesian for overlap check
  const x1 = box1.radius * Math.cos(box1.angle);
  const y1 = box1.radius * Math.sin(box1.angle);
  const x2 = box2.radius * Math.cos(box2.angle);
  const y2 = box2.radius * Math.sin(box2.angle);
  
  const dx = x2 - x1;
  const dy = y2 - y1;
  const dist = Math.sqrt(dx * dx + dy * dy);
  
  const minDist = (box1.width + box2.width) / 2 + minSeparation;
  return dist < minDist;
}

/**
 * Group residues by type for outer ring organization
 */
export function groupResiduesByType(residues) {
  const groups = {
    charged_negative: [],
    charged_positive: [],
    polar: [],
    hydrophobic: [],
    glycine: [],
    nucleic: [],
    other: [],
  };
  
  residues.forEach(res => {
    const cls = res.class || "other";
    if (cls === "negative") {
      groups.charged_negative.push(res);
    } else if (cls === "positive") {
      groups.charged_positive.push(res);
    } else if (cls === "polar") {
      groups.polar.push(res);
    } else if (cls === "hydrophobic") {
      groups.hydrophobic.push(res);
    } else if (cls === "glycine") {
      groups.glycine.push(res);
    } else if (cls === "nucleic") {
      groups.nucleic.push(res);
    } else {
      groups.other.push(res);
    }
  });
  
  return groups;
}

/**
 * Calculate angular spacing for a group of residues
 */
export function calculateAngularSpacing(residues, minAngle = 0.15) {
  if (residues.length === 0) return [];
  
  // Sort by angle
  const sorted = [...residues].sort((a, b) => (a.angle || 0) - (b.angle || 0));
  
  // Calculate required spacing
  const totalAngle = 2 * Math.PI;
  const requiredSpacing = Math.max(
    minAngle,
    totalAngle / sorted.length * 1.2 // 20% extra spacing
  );
  
  return sorted.map((res, idx) => ({
    ...res,
    targetAngle: idx * requiredSpacing,
  }));
}

/**
 * Apply angular clustering for residues within 1 Å distance
 */
export function applyAngularClustering(residues, distanceThreshold = 1.0) {
  // Group residues that are close in 3D space
  const clusters = [];
  const processed = new Set();
  
  residues.forEach((res, idx) => {
    if (processed.has(idx)) return;
    
    const cluster = [res];
    processed.add(idx);
    
    residues.forEach((other, otherIdx) => {
      if (processed.has(otherIdx)) return;
      
      // Check if residues are close in 3D distance
      const dist = res.dist || 0;
      const otherDist = other.dist || 0;
      const distDiff = Math.abs(dist - otherDist);
      
      if (distDiff < distanceThreshold) {
        cluster.push(other);
        processed.add(otherIdx);
      }
    });
    
    if (cluster.length > 0) {
      clusters.push(cluster);
    }
  });
  
  return clusters;
}

/**
 * Radial push-out collision resolution
 */
export function resolveCollisions(residues, ligandCenter, minRadius = 200, maxRadius = 400) {
  const resolved = residues.map(r => ({ ...r }));
  const maxIterations = 50;
  const pushOutFactor = 1.1;
  const minAngularSeparation = 0.12; // ~7 degrees
  
  for (let iter = 0; iter < maxIterations; iter++) {
    let hasOverlap = false;
    
    // Sort by angle
    resolved.sort((a, b) => (a.angle || 0) - (b.angle || 0));
    
    for (let i = 0; i < resolved.length; i++) {
      const res1 = resolved[i];
      const box1 = getResidueBoundingBox(res1, res1.radius || minRadius);
      
      for (let j = i + 1; j < resolved.length; j++) {
        const res2 = resolved[j];
        const box2 = getResidueBoundingBox(res2, res2.radius || minRadius);
        
        // Check angular separation first
        const angleDiff = Math.abs((res2.angle || 0) - (res1.angle || 0));
        const normalizedAngleDiff = Math.min(angleDiff, 2 * Math.PI - angleDiff);
        
        if (normalizedAngleDiff < minAngularSeparation) {
          // Too close angularly - push one out radially
          if (res1.radius < res2.radius) {
            res1.radius = Math.min((res1.radius || minRadius) * pushOutFactor, maxRadius);
          } else {
            res2.radius = Math.min((res2.radius || minRadius) * pushOutFactor, maxRadius);
          }
          hasOverlap = true;
        } else if (boxesOverlap(box1, box2)) {
          // Bounding boxes overlap - push out radially
          if (res1.radius < res2.radius) {
            res1.radius = Math.min((res1.radius || minRadius) * pushOutFactor, maxRadius);
          } else {
            res2.radius = Math.min((res2.radius || minRadius) * pushOutFactor, maxRadius);
          }
          hasOverlap = true;
        }
      }
    }
    
    if (!hasOverlap) break;
  }
  
  return resolved;
}

/**
 * Compute interaction anchor point on ligand for a residue
 * Uses the closest ligand atom from interactions, or estimates from residue angle
 */
export function computeInteractionAnchor(residue, interactions, ligandAtoms, ligandCenter) {
  // Find interactions for this residue
  const residueInteractions = interactions.filter(
    it => it.residue === `${residue.resname}${residue.resid}`
  );
  
  if (residueInteractions.length > 0 && ligandAtoms) {
    // Use the first interaction's ligand atom as anchor
    const firstInteraction = residueInteractions[0];
    const ligandAtomIdx = firstInteraction.ligand_atom_index;
    if (ligandAtoms[ligandAtomIdx]) {
      return {
        x: ligandAtoms[ligandAtomIdx].x,
        y: ligandAtoms[ligandAtomIdx].y,
        angle: Math.atan2(
          ligandAtoms[ligandAtomIdx].y - ligandCenter[1],
          ligandAtoms[ligandAtomIdx].x - ligandCenter[0]
        ),
      };
    }
  }
  
  // Fallback: estimate anchor from residue's biological angle
  const angle = residue.angle || 0;
  const estimatedRadius = 80; // Approximate ligand radius
  return {
    x: ligandCenter[0] + estimatedRadius * Math.cos(angle),
    y: ligandCenter[1] + estimatedRadius * Math.sin(angle),
    angle: angle,
  };
}

/**
 * Cluster residues by angular region (sector clustering)
 */
export function clusterResiduesBySector(residues, interactions, ligandAtoms, ligandCenter, numSectors = 8) {
  // Compute interaction anchors for all residues
  const residuesWithAnchors = residues.map(res => {
    const anchor = computeInteractionAnchor(res, interactions, ligandAtoms, ligandCenter);
    return {
      ...res,
      anchorAngle: anchor.angle,
      anchorX: anchor.x,
      anchorY: anchor.y,
    };
  });
  
  // Normalize angles to [0, 2π]
  residuesWithAnchors.forEach(res => {
    res.anchorAngle = ((res.anchorAngle % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI);
  });
  
  // Sort by anchor angle
  residuesWithAnchors.sort((a, b) => a.anchorAngle - b.anchorAngle);
  
  // Cluster into sectors
  const sectorSize = (2 * Math.PI) / numSectors;
  const clusters = [];
  let currentCluster = [];
  let currentSector = 0;
  
  residuesWithAnchors.forEach(res => {
    const sector = Math.floor(res.anchorAngle / sectorSize);
    
    if (sector === currentSector) {
      currentCluster.push(res);
    } else {
      if (currentCluster.length > 0) {
        clusters.push({
          sector: currentSector,
          residues: [...currentCluster],
          minAngle: currentCluster[0].anchorAngle,
          maxAngle: currentCluster[currentCluster.length - 1].anchorAngle,
        });
      }
      currentCluster = [res];
      currentSector = sector;
    }
  });
  
  // Don't forget the last cluster
  if (currentCluster.length > 0) {
    clusters.push({
      sector: currentSector,
      residues: [...currentCluster],
      minAngle: currentCluster[0].anchorAngle,
      maxAngle: currentCluster[currentCluster.length - 1].anchorAngle,
    });
  }
  
  return clusters;
}

/**
 * Allocate arc ranges per cluster with gaps between clusters
 */
export function allocateArcRanges(clusters, minGap = 0.1) {
  if (clusters.length === 0) return [];
  
  const totalAngle = 2 * Math.PI;
  const totalGaps = (clusters.length - 1) * minGap;
  const availableAngle = totalAngle - totalGaps;
  
  // Calculate total "weight" (number of residues) across all clusters
  const totalResidues = clusters.reduce((sum, cluster) => sum + cluster.residues.length, 0);
  
  // Allocate arcs proportionally to cluster size
  let currentAngle = 0;
  const allocatedClusters = clusters.map((cluster, idx) => {
    const clusterWeight = cluster.residues.length / totalResidues;
    const arcSize = availableAngle * clusterWeight;
    
    const arcStart = currentAngle;
    const arcEnd = currentAngle + arcSize;
    
    currentAngle = arcEnd + (idx < clusters.length - 1 ? minGap : 0);
    
    return {
      ...cluster,
      arcStart: arcStart,
      arcEnd: arcEnd,
      arcSize: arcSize,
    };
  });
  
  return allocatedClusters;
}

/**
 * Sector-based constrained pocket layout - Maestro/LigPlot style
 */
export function computeRadialLayout(residues, ligandCenter, interactions = [], ligandAtoms = null, options = {}) {
  const {
    minRadius = 200,
    maxRadius = 400,
    minAngularSpacing = 0.12, // ~6.9 degrees minimum within arc
    numSectors = 8,
    interClusterGap = 0.15, // Gap between clusters
  } = options;
  
  if (!residues || residues.length === 0) {
    return [];
  }
  
  // Step 1: Cluster residues by angular region (sector clustering)
  const clusters = clusterResiduesBySector(
    residues,
    interactions,
    ligandAtoms,
    ligandCenter,
    numSectors
  );
  
  // Step 2: Allocate arc ranges per cluster with gaps
  const allocatedClusters = allocateArcRanges(clusters, interClusterGap);
  
  // Step 3: Place residues within their assigned arcs
  const layoutResidues = [];
  
  allocatedClusters.forEach(cluster => {
    const clusterResidues = cluster.residues;
    const arcSize = cluster.arcSize;
    const arcStart = cluster.arcStart;
    
    if (clusterResidues.length === 0) return;
    
    // Evenly distribute residues within the arc
    const spacing = Math.max(
      arcSize / clusterResidues.length,
      minAngularSpacing
    );
    
    clusterResidues.forEach((res, idx) => {
      // Map 3D distance to 2D radius (variable based on pocket contour)
      const dist3D = res.dist || 0;
      const normalizedDist = Math.min(Math.max(dist3D, 0), 20);
      const baseRadius = minRadius + (normalizedDist / 20) * (maxRadius - minRadius);
      
      // Assign angle within the arc
      const assignedAngle = arcStart + idx * spacing;
      
      layoutResidues.push({
        ...res,
        angle: assignedAngle,
        radius: baseRadius,
        clusterSector: cluster.sector,
        anchorAngle: res.anchorAngle,
        anchorX: res.anchorX,
        anchorY: res.anchorY,
      });
    });
  });
  
  // Step 4: Resolve collisions with iterative relaxation
  const resolvedResidues = resolveCollisions(layoutResidues, ligandCenter, minRadius, maxRadius);
  
  // Step 5: Convert polar to cartesian coordinates
  const finalResidues = resolvedResidues.map(res => {
    const x = ligandCenter[0] + res.radius * Math.cos(res.angle);
    const y = ligandCenter[1] + res.radius * Math.sin(res.angle);
    
    return {
      ...res,
      x: x,
      y: y,
    };
  });
  
  return finalResidues;
}

/**
 * Generate quadratic Bézier path between two points
 */
export function createQuadraticBezierPath(x1, y1, x2, y2, curvature = 0.3) {
  const dx = x2 - x1;
  const dy = y2 - y1;
  const dist = Math.sqrt(dx * dx + dy * dy);
  
  // Control point perpendicular to the line (midpoint with offset)
  const midX = (x1 + x2) / 2;
  const midY = (y1 + y2) / 2;
  const angle = Math.atan2(dy, dx);
  const perpAngle = angle + Math.PI / 2;
  const offset = dist * curvature;
  
  const cpx = midX + Math.cos(perpAngle) * offset;
  const cpy = midY + Math.sin(perpAngle) * offset;
  
  return `M ${x1} ${y1} Q ${cpx} ${cpy}, ${x2} ${y2}`;
}

/**
 * Generate two-segment edge routing: residue → sector anchor curve → ligand atom
 * Uses quadratic Bézier curves for smooth routing
 * Sector anchor is placed based on the residue's interaction anchor point
 */
export function createTwoSegmentPath(
  residueX, 
  residueY, 
  ligandAtomX, 
  ligandAtomY, 
  ligandCenter,
  interactionAnchorX = null,
  interactionAnchorY = null,
  anchorRadius = 120
) {
  // Use interaction anchor if provided, otherwise calculate from residue angle
  let sectorAnchorX, sectorAnchorY;
  
  if (interactionAnchorX !== null && interactionAnchorY !== null) {
    // Use the interaction anchor point as the sector anchor
    const angleToAnchor = Math.atan2(
      interactionAnchorY - ligandCenter[1],
      interactionAnchorX - ligandCenter[0]
    );
    sectorAnchorX = ligandCenter[0] + anchorRadius * Math.cos(angleToAnchor);
    sectorAnchorY = ligandCenter[1] + anchorRadius * Math.sin(angleToAnchor);
  } else {
    // Fallback: calculate from residue position
    const angleToResidue = Math.atan2(residueY - ligandCenter[1], residueX - ligandCenter[0]);
    sectorAnchorX = ligandCenter[0] + anchorRadius * Math.cos(angleToResidue);
    sectorAnchorY = ligandCenter[1] + anchorRadius * Math.sin(angleToResidue);
  }
  
  // First segment: residue → sector anchor (quadratic Bézier)
  const seg1 = createQuadraticBezierPath(residueX, residueY, sectorAnchorX, sectorAnchorY, 0.25);
  
  // Second segment: sector anchor → ligand atom (quadratic Bézier)
  const seg2 = createQuadraticBezierPath(sectorAnchorX, sectorAnchorY, ligandAtomX, ligandAtomY, 0.25);
  
  return {
    path: `${seg1} ${seg2.substring(1)}`, // Remove duplicate M from second segment
    anchorX: sectorAnchorX,
    anchorY: sectorAnchorY,
  };
}

/**
 * Generate circular path for consecutive residue group
 */
export function createCircularPath(residues, ligandCenter, isInnerCircle = false) {
  if (residues.length < 2) return null;
  
  const pathPoints = residues.map(r => ({
    x: r.x,
    y: r.y,
    angle: Math.atan2(r.y - ligandCenter[1], r.x - ligandCenter[0]),
  }));
  
  // Sort by angle for smooth circular flow
  pathPoints.sort((a, b) => a.angle - b.angle);
  
  const curveDirection = isInnerCircle ? -1 : 1;
  const baseOffset = isInnerCircle ? 25 : 35;
  
  let pathData = `M ${pathPoints[0].x} ${pathPoints[0].y}`;
  
  for (let i = 0; i < pathPoints.length; i++) {
    const current = pathPoints[i];
    const next = pathPoints[(i + 1) % pathPoints.length];
    
    const dx = next.x - current.x;
    const dy = next.y - current.y;
    const dist = Math.sqrt(dx * dx + dy * dy);
    const offset = Math.min(dist * 0.4, baseOffset) * curveDirection;
    
    const angle = Math.atan2(dy, dx);
    const perpAngle = angle + Math.PI / 2;
    
    const cp1x = current.x + Math.cos(perpAngle) * offset;
    const cp1y = current.y + Math.sin(perpAngle) * offset;
    const cp2x = next.x + Math.cos(perpAngle) * offset;
    const cp2y = next.y + Math.sin(perpAngle) * offset;
    
    if (i === pathPoints.length - 1 && pathPoints.length > 2) {
      // Close the circle
      const first = pathPoints[0];
      const dxClose = first.x - current.x;
      const dyClose = first.y - current.y;
      const angleClose = Math.atan2(dyClose, dxClose);
      const perpAngleClose = angleClose + Math.PI / 2;
      const offsetClose = Math.min(
        Math.sqrt(dxClose * dxClose + dyClose * dyClose) * 0.4,
        baseOffset
      ) * curveDirection;
      
      const cp1CloseX = current.x + Math.cos(perpAngleClose) * offsetClose;
      const cp1CloseY = current.y + Math.sin(perpAngleClose) * offsetClose;
      const cp2CloseX = first.x + Math.cos(perpAngleClose) * offsetClose;
      const cp2CloseY = first.y + Math.sin(perpAngleClose) * offsetClose;
      
      pathData += ` C ${cp1CloseX} ${cp1CloseY}, ${cp2CloseX} ${cp2CloseY}, ${first.x} ${first.y} Z`;
    } else {
      pathData += ` C ${cp1x} ${cp1y}, ${cp2x} ${cp2y}, ${next.x} ${next.y}`;
    }
  }
  
  return pathData;
}
