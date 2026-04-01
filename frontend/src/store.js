import { configureStore } from "@reduxjs/toolkit";
import diagramReducer from "./features/diagramSlice";

const store = configureStore({
  reducer: {
    diagram: diagramReducer,
  },
});

export default store;