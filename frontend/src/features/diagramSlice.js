import { createSlice, createAsyncThunk } from "@reduxjs/toolkit";
import api from "../api";

export const generateDiagram = createAsyncThunk(
  "diagram/generate",
  async ({ pdbFile, sdfFile, ligandResname, pocketRadius, modelIndex }) => {
    const form = new FormData();
    form.append("protein_ligand_pdb", pdbFile);
    // Only append SDF if provided (optional)
    if (sdfFile) {
      form.append("ligand_sdf", sdfFile);
    }

    const params = {
      pocket_radius: pocketRadius,
      model_index: modelIndex || 1,
    };
    // Only add ligand_resname if provided (not empty)
    if (ligandResname && ligandResname.trim() !== "") {
      params.ligand_name = ligandResname.trim();
    }

    const response = await api.post("/lid/generate", form, {
      params: params,
      headers: { "Content-Type": "multipart/form-data" },
    });

    return response.data;
  }
);

const diagramSlice = createSlice({
  name: "diagram",
  initialState: {
    loading: false,
    data: null,
    error: null,
  },
  reducers: {},
  extraReducers: (builder) => {
    builder
      .addCase(generateDiagram.pending, (state) => {
        state.loading = true;
        state.error = null;
      })
      .addCase(generateDiagram.fulfilled, (state, action) => {
        state.loading = false;
        state.data = action.payload;
      })
      .addCase(generateDiagram.rejected, (state, action) => {
        state.loading = false;
        state.error = action.error.message;
      });
  },
});

export default diagramSlice.reducer;