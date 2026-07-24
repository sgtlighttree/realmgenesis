import React from 'react';
import ReactDOM from 'react-dom/client';
import App from './App';
import DesignShell from './components/shell/DesignShell';
import './index.css';

const rootElement = document.getElementById('root');
if (!rootElement) {
  throw new Error("Could not find root element to mount to");
}

// ?shell=1 mounts the F1 layout prototype instead of the live app.
const useShell = new URLSearchParams(window.location.search).has('shell');

const root = ReactDOM.createRoot(rootElement);
root.render(
  <React.StrictMode>
    {useShell ? <DesignShell /> : <App />}
  </React.StrictMode>
);