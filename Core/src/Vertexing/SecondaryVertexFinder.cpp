// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Vertexing/SecondaryVertexFinder.hpp"

auto Acts::SecondaryVertexFinder::find(
    const std::vector<InputTrack>& trackVector,
    const VertexingOptions& vertexingOptions,
    IVertexFinder::State& anyState) const -> Result<std::vector<Vertex>> {
  auto& state = anyState.as<State>();
  // Original tracks
  const std::vector<InputTrack>& origTracks = trackVector;
  // Tracks for seeding
  std::vector<InputTrack> seedTracks = trackVector;

  // List of vertices to be filled below
  std::vector<Vertex> vertexCollection;

  while (seedTracks.size() > 1) {
    auto seedRes = getVertexSeed(state, seedTracks, vertexingOptions);

    if (!seedRes.ok()) {
      return seedRes.error();
    }
    const auto& seedOptional = *seedRes;

    if (!seedOptional.has_value()) {
      ACTS_DEBUG("No more seed found. Break and stop secondary vertex finding.");
      break;
    }
    const auto& seedVertex = *seedOptional;

    /// End seeding
    /// Now take only tracks compatible with current seed
    // Tracks used for the fit in this iteration
    std::vector<InputTrack> tracksToFit;
    std::vector<InputTrack> tracksToFitSplitVertex;

    // Fill vector with tracks to fit, only compatible with seed:
    auto res = fillTracksToFit(seedTracks, seedVertex, tracksToFit,
                               tracksToFitSplitVertex, vertexingOptions, state);

    if (!res.ok()) {
      return res.error();
    }

    ACTS_DEBUG("Number of tracks used for fit: " << tracksToFit.size());

    /// Begin vertex fit
    Vertex currentVertex;
    Vertex currentSplitVertex;

    if (vertexingOptions.useConstraintInFit && !tracksToFit.empty()) {
      auto fitResult = m_cfg.vertexFitter.fit(tracksToFit, vertexingOptions,
                                              state.fieldCache);
      if (fitResult.ok()) {
        currentVertex = std::move(*fitResult);
      } else {
        return fitResult.error();
      }
    } else if (!vertexingOptions.useConstraintInFit && tracksToFit.size() > 1) {
      auto fitResult = m_cfg.vertexFitter.fit(tracksToFit, vertexingOptions,
                                              state.fieldCache);
      if (fitResult.ok()) {
        currentVertex = std::move(*fitResult);
      } else {
        return fitResult.error();
      }
    }
    if (m_cfg.createSplitVertices && tracksToFitSplitVertex.size() > 1) {
      auto fitResult = m_cfg.vertexFitter.fit(
          tracksToFitSplitVertex, vertexingOptions, state.fieldCache);
      if (fitResult.ok()) {
        currentSplitVertex = std::move(*fitResult);
      } else {
        return fitResult.error();
      }
    }
    /// End vertex fit
    ACTS_DEBUG("Vertex position after fit: "
               << currentVertex.fullPosition().transpose());

    // Number degrees of freedom
    double ndf = currentVertex.fitQuality().second;
    double ndfSplitVertex = currentSplitVertex.fitQuality().second;

    // Number of significant tracks
    int nTracksAtVertex = countSignificantTracks(currentVertex);
    int nTracksAtSplitVertex = countSignificantTracks(currentSplitVertex);

    bool isGoodVertex = ((!vertexingOptions.useConstraintInFit && ndf > 0 &&
                          nTracksAtVertex >= 2) ||
                         (vertexingOptions.useConstraintInFit && ndf > 3 &&
                          nTracksAtVertex >= 2));

    if (!isGoodVertex) {
      removeTracks(tracksToFit, seedTracks);
    } else {
      if (m_cfg.reassignTracksAfterFirstFit && (!m_cfg.createSplitVertices)) {
        // vertex is good vertex here
        // but add tracks which may have been missed

        auto result = reassignTracksToNewVertex(
            vertexCollection, currentVertex, tracksToFit, seedTracks,
            origTracks, vertexingOptions, state);
        if (!result.ok()) {
          return result.error();
        }
        isGoodVertex = *result;

      }  // end reassignTracksAfterFirstFit case
         // still good vertex? might have changed in the meanwhile
      if (isGoodVertex) {
        removeUsedCompatibleTracks(currentVertex, tracksToFit, seedTracks,
                                   vertexingOptions, state);

        ACTS_DEBUG(
            "Number of seed tracks after removal of compatible tracks "
            "and outliers: "
            << seedTracks.size());
      }
    }  // end case if good vertex

    // now splitvertex
    bool isGoodSplitVertex = false;
    if (m_cfg.createSplitVertices) {
      isGoodSplitVertex = (ndfSplitVertex > 0 && nTracksAtSplitVertex >= 2);

      if (!isGoodSplitVertex) {
        removeTracks(tracksToFitSplitVertex, seedTracks);
      } else {
        removeUsedCompatibleTracks(currentSplitVertex, tracksToFitSplitVertex,
                                   seedTracks, vertexingOptions, state);
      }
    }
    // Now fill vertex collection with vertex
    if (isGoodVertex) {
      vertexCollection.push_back(currentVertex);
    }
    if (isGoodSplitVertex && m_cfg.createSplitVertices) {
      vertexCollection.push_back(currentSplitVertex);
    }
  }

  return vertexCollection;
}
