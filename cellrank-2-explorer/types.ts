import * as THREE from 'three';

export type AppMode = 'MENU' | 'EXPLORER' | 'FORMULAS';

export enum GameStage {
  INTRO = 0,
  KNN_GRAPH = 1,
  KERNEL_BIAS = 2,
  TRANSITION_MATRIX = 3,
  RANDOM_WALK = 4,
  MACROSTATES = 5,
}

export type CellData = {
  id: number;
  position: THREE.Vector3;
  color: string;
  pseudotime: number; // 0 to 1
  potency: number; // 0 to 1 (CytoTRACE score)
  velocity: THREE.Vector3; // For velocity kernel
  neighbors: number[]; // IDs of neighbors
  type: 'Stem' | 'Progenitor' | 'TypeA' | 'TypeB';
};

export type KernelType = 'Pseudotime' | 'Velocity' | 'CytoTRACE' | 'RealTime' | 'Combined';

export interface SimulationState {
  cells: CellData[];
  activeCell: number | null;
  walkPath: number[];
  terminalStates: number[]; // IDs of representative terminal cells
}