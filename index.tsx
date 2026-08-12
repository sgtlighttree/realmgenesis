import React from 'react';
import ReactDOM from 'react-dom/client';
import App from './App';
import ShellApp from './components/shell/ShellApp';
import DesignShell from './components/shell/DesignShell';
import './index.css';

const rootElement = document.getElementById('root');
if (!rootElement) {
  throw new Error("Could not find root element to mount to");
}

// Entry routing (F1 shell is now the default after the parity smoke passed):
//   ?shell=stub    → DesignShell (F1 layout prototype, stub panels)
//   ?shell=classic → classic App (legacy fork, kept until fully retired)
//   otherwise      → ShellApp (F1 redesign, real data) — incl. old ?shell=1 links
const shellParam = new URLSearchParams(window.location.search).get('shell');

const Root =
  shellParam === 'stub' ? <DesignShell />
  : shellParam === 'classic' ? <App />
  : <ShellApp />;

const root = ReactDOM.createRoot(rootElement);
root.render(
  <React.StrictMode>
    {Root}
  </React.StrictMode>
);