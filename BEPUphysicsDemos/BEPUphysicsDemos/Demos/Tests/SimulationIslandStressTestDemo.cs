﻿using System;
using BEPUphysics.Collidables;
using BEPUphysics.Entities.Prefabs;
using BEPUphysics.MathExtensions;
using Microsoft.Xna.Framework;
using BEPUphysics.Settings;
using BEPUphysics.CollisionTests.CollisionAlgorithms;
using BEPUphysics.CollisionRuleManagement;
using BEPUphysics.NarrowPhaseSystems.Pairs;
using BEPUphysics.CollisionTests.Manifolds;
using System.Diagnostics;
using BEPUphysics.BroadPhaseSystems.SortAndSweep;
using BEPUphysics.NarrowPhaseSystems;

namespace BEPUphysicsDemos.Demos.Tests
{
    /// <summary>
    /// Spheres fall onto a large terrain.  Try driving around on it!
    /// </summary>
    public class SimulationIslandStressTestDemo : StandardDemo
    {
        /// <summary>
        /// Constructs a new demo.
        /// </summary>
        /// <param name="game">Game owning this demo.</param>
        public SimulationIslandStressTestDemo(DemosGame game)
            : base(game)
        {

            Space.Add(new Box(new Vector3(0, 0, 0), 500, 10, 500));
            //MotionSettings.DefaultPositionUpdateMode = BEPUphysics.PositionUpdating.PositionUpdateMode.Continuous;

            NarrowPhaseHelper.Factories.SphereSphere.EnsureCount(10000);
            NarrowPhaseHelper.Factories.BoxSphere.EnsureCount(3000);
            //Space.BroadPhase = new Grid2DSortAndSweep(Space.ThreadManager);

            //ConfigurationHelper.ApplySuperSpeedySettings(Space);
            //Space.ForceUpdater.Gravity = new Vector3();
            Space.Solver.IterationLimit = 1;
            //for (int i = 0; i < 20; i++)
            //{
            //    for (int j = 0; j < 10; j++)
            //    {
            //        for (int k = 0; k < 20; k++)
            //        {
            //            Space.Add(new Sphere(new Vector3(0 + i * 3, 20 + j * 3, 0 + k * 3), 0.5f, 1)
            //            {
            //                //LocalInertiaTensorInverse = new Matrix3X3()
            //                //Tag = "noDisplayObject",
            //                IsAlwaysActive = true
            //            });
            //        }
            //    }
            //}


            Random rand = new Random();


            float width = 30;
            float height = 200;
            float length = 30; 
            for (int i = 0; i < 3000; i++)
            {
                Vector3 position =
                    new Vector3((float)rand.NextDouble() * width - width * .5f,
                        (float)rand.NextDouble() * height + 20,
                        (float)rand.NextDouble() * length - length * .5f);
                Space.Add(new Sphere(position, 1, 1) { Tag = "noDisplayObject" }
                );
            }





            game.Camera.Position = new Vector3(0, 30, 20);

            ////Pre-simulate.
            //for (int i = 0; i < 150; i++)
            //{
            //    Space.Update();
            //}

            //int numRuns = 10;
            ////Space.BeforeNarrowPhaseUpdateables.Enabled = false;
            ////Space.DuringForcesUpdateables.Enabled = false;
            ////Space.EndOfTimeStepUpdateables.Enabled = false;
            ////Space.EndOfFrameUpdateables.Enabled = false;

            ////Space.EntityStateWriteBuffer.Enabled = false;
            ////Space.BufferedStates.Enabled = false;
            ////Space.DeactivationManager.Enabled = false;
            ////Space.SpaceObjectBuffer.Enabled = false;
            ////Space.DeferredEventDispatcher.Enabled = false;
            ////Space.BoundingBoxUpdater.Enabled = false;

            ////Space.BroadPhase.Enabled = false;
            ////Space.NarrowPhase.Enabled = false;
            ////Space.ForceUpdater.Enabled = false;
            //////Space.Solver.Enabled = false;
            ////Space.PositionUpdater.Enabled = false;
            //for (int i = 0; i < numRuns; i++)
            //{
            //    Space.Update();
            //}

        }




        /// <summary>
        /// Gets the name of the simulation.
        /// </summary>
        public override string Name
        {
            get { return "Simulation Island Stress Test"; }
        }
    }
}